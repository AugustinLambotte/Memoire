CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   N_CALIB       	N_HISTORY             
   title         Argo float vertical profile    institution       AOML   source        
Argo float     history       2018-10-06T00:40:27Z creation      
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
_FillValue                    �tArgo profile    3.1 1.2 19500101000000  20181006004027  20190405134400  5904772 US ARGO PROJECT                                                 STEPHEN RISER                                                   PRES            TEMP            PSAL               <A   AO  6627                            2C  D   APEX                            7392                            062915                          846 @�@���1   @�@���Z~@N������@��^5?}1   GPS     Primary sampling: mixed [deeper than nominal 985dbar: discrete; nominal 985dbar to surface: 2dbar-bin averaged]                                                                                                                                                    <A   B   B   @9��@�  @�  A   A   A@  A`  A�  A�  A�  A�33A�  A�33A�ffA�  B ffBffB��B  B   B(  B0  B8  B@  BH  BO��BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� DmfDm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Dt�fDy��D�fD�<�D���D�ٚD��D�33D���D�� D�	�D��3D��fDǼ�D�	�D�0 Dڙ�D๚D��3D�,�D�|�D��31111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @q�@�(�@�(�A{A.{AN{An{A�
=A�
=A�
=A�=pA�
=A�=pA�p�A�
=B�B�B�B�B#�B+�B3�B;�BC�BK�BS�B[�Bc�Bk�Bs�B{�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�C �HC�HC�HC�HC�HC
�HC�HC�HC�HC�HC�HC�HC�HC�HC�HC�HC �HC"�HC$�HC&�HC(�HC*�HC,�HC.�HC0�HC2�HC4�HC6�HC8�HC:�HC<�HC>�HC@�HCB�HCD�HCF�HCH�HCJ�HCL�HCN�HCP�HCR�HCT�HCV�HCX�HCZ�HC\�HC^�HC`�HCb�HCd�HCf�HCh�HCj�HCl�HCn�HCp�HCr�HCt�HCv�HCx�HCz�HC|�HC~�HC�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�p�C�c�D 8RD �RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD	8RD	�RD
8RD
�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD8RD�RD 8RD �RD!8RD!�RD"8RD"�RD#8RD#�RD$8RD$�RD%8RD%�RD&8RD&�RD'8RD'�RD(8RD(�RD)8RD)�RD*8RD*�RD+8RD+�RD,8RD,�RD-8RD-�RD.8RD.�RD/8RD/�RD08RD0�RD18RD1�RD28RD2�RD38RD3�RD48RD4�RD58RD5�RD68RD6�RD78RD7�RD88RD8�RD98RD9�RD:8RD:�RD;8RD;�RD<8RD<�RD=8RD=�RD>8RD>�RD?8RD?�RD@8RD@�RDA8RDA�RDB8RDB�RDC8RDC�RDD8RDD�RDE8RDE�RDF8RDF�RDG8RDG�RDH8RDH�RDI8RDI�RDJ8RDJ�RDK8RDK�RDL8RDL�RDM8RDM�RDN8RDN�RDO8RDO�RDP8RDP�RDQ8RDQ�RDR8RDR�RDS8RDS�RDT8RDT�RDU8RDU�RDV8RDV�RDW8RDW�RDX8RDX�RDY8RDY�RDZ8RDZ�RD[8RD[�RD\8RD\�RD]8RD]�RD^8RD^�RD_8RD_�RD`8RD`�RDa8RDa�RDb8RDb�RDc8RDc�RDd8RDd�RDe8RDe�RDf8RDf�RDg8RDg�RDh8RDh�RDi8RDi�RDj8RDj�RDk8RDk�RDl8RDl�RDm>�Dm�RDn8RDn�RDo8RDo�RDp8RDp�RDq8RDq�RDr8RDr�RDs8RDs�RDt8RDt�RDu�Dy��D�"�D�X�D���D���D�(�D�O\D���D��)D�%�D��\D���D���D�%�D�L)Dڵ�D���D�\D�H�D��D��\1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���@���@�
=@�
=@�
=@�
=@�
=@�@�
=@�
=@�o@�o@�o@�o@�o@�"�@��@��@�+@�"�@�+@�+@�+@�+@�+@�+@�+@�+@�+@�+@�+@�+@�+@�+@�+@�+@�33@�;d@�33@�C�@�C�@�K�@�C�@�C�@�;d@�
=@���@��@��y@�ȴ@���@���@��+@�n�@�^5@�E�@�$�@�@���@���@���@�hs@�O�@��@�V@�V@�V@���@��`@��j@��u@�z�@�r�@�j@�r�@�j@�j@�bN@�bN@�bN@�j@�j@�j@�r�@�r�@�r�@��@��u@���@��@��9@��j@���@���@���@���@���@���@���@��@��#@��@�/@���@�9X@��;@��w@���@�t�@�\)@�dZ@�C�@�+@�+@�"�@�@���@�ff@�M�@�=q@�5?@�=q@�5?@�5?@�-@�-@��@�J@�@�x�@�`B@�?}@�&�@��@���@��u@�9X@�(�@��m@��
@��w@���@�;d@��@���@���@�ff@�5?@�{@��@��^@�&�@���@��@��u@��u@�z�@�bN@�1'@��@|�@;d@+@�@~��@~ȴ@~ȴ@~��@~v�@~V@~$�@}�@}�T@}��@}��@}`B@|��@|��@|I�@{��@{33@{@z�H@z��@z�\@zn�@zM�@zM�@zM�@z�@y��@y�7@yG�@y%@x�`@x�`@x�`@x�`@x�`@x��@x�@xQ�@xQ�@x �@w�@w�@w�@w�;@w��@w��@w�w@wl�@wK�@v��@v�@v�R@v��@v�+@vv�@v@u�-@u�h@u�h@u�h@u�@up�@uO�@u/@t��@t�/@tz�@s��@sC�@so@s@r�H@r��@r��@r�!@r�!@r�!@r�!@r�!@r�!@r�\@r-@q��@q��@q�^@q�^@q�^@q�7@qhs@qG�@qX@qX@qhs@qX@q%@pĜ@p�u@pbN@pQ�@p1'@pb@o�@o�@o�@o�@o\)@o\)@o\)@o\)@o\)@o\)@oK�@o�@n��@nȴ@nȴ@n�R@n��@n��@n�+@n�+@nV@nV@n$�@m�T@m�@mp�@mO�@m/@m�@l�@l�/@l��@l�j@l�j@l�@l��@l�D@l�D@lj@l1@k��@kƨ@k�F@k�F@k�F@k��@k�@kt�@kdZ@kC�@k"�@ko@j�@j��@j�!@j��@j��@j��@j��@j�\@j~�@jn�@jn�@jn�@j^5@j^5@j^5@jM�@j=q@j-@i��@i��@i�@i�#@i��@i��@i��@ix�@ihs@ihs@iX@iX@iX@iG�@i7L@i7L@i&�@h��@h�`@h�`@h�`@h�`@h�`@h�`@h��@h��@hĜ@hĜ@hĜ@h��@hĜ@h�9@h�9@h�9@h��@h�u@h�u@hr�@hbN@hbN@hbN@hbN@hbN@hr�@hr�@hr�@hr�@h�@h�@h�@h�@h�@h�`@h��@h�`@h�`@hĜ@h�9@h�u@hr�@hbN@hQ�@hQ�@hQ�@hQ�@hA�@hA�@hA�@h1'@h1'@h �@hb@hb@h �@h �@h �@h �@hb@hb@h �@hb@hb@hb@h  @hb@h �@h �@hb@h  @g�@g�;@g�;@g��@g�w@g�w@g�@g�@g��@g��@g�P@g�P@g|�@gl�@gK�@gK�@g;d@g;d@g+@g;d@g;d@g+@g�@g
=@f��@f��@f�y@f�y@fȴ@f�R@f��@f��@f�+@f�+@fv�@fff@fV@fE�@f5?@f$�@f{@f@e�T@e��@e@e�-@e�-@e��@e�h@e�@ep�@ep�@e`B@eO�@e?}@e?}@e?}@e/@e/@e/@e�@eV@eV@eV@d��@d��@d�@d9X@c"�@a��@_�;@_+@_
=@^ff@_�P@`  @co@d(�@f��@h  @g\)@f@ct�@a��@`��@]�T@Z~�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 @���@���@�
=@�
=@�
=@�
=@�
=@�@�
=@�
=@�o@�o@�o@�o@�o@�"�@��@��@�+@�"�@�+@�+@�+@�+@�+@�+@�+@�+@�+@�+@�+@�+@�+@�+@�+@�+@�33@�;d@�33@�C�@�C�@�K�@�C�@�C�@�;d@�
=@���@��@��y@�ȴ@���@���@��+@�n�@�^5@�E�@�$�@�@���@���@���@�hs@�O�@��@�V@�V@�V@���@��`@��j@��u@�z�@�r�@�j@�r�@�j@�j@�bN@�bN@�bN@�j@�j@�j@�r�@�r�@�r�@��@��u@���@��@��9@��j@���@���@���@���@���@���@���@��@��#@��@�/@���@�9X@��;@��w@���@�t�@�\)@�dZ@�C�@�+@�+@�"�@�@���@�ff@�M�@�=q@�5?@�=q@�5?@�5?@�-@�-@��@�J@�@�x�@�`B@�?}@�&�@��@���@��u@�9X@�(�@��m@��
@��w@���@�;d@��@���@���@�ff@�5?@�{@��@��^@�&�@���@��@��u@��u@�z�@�bN@�1'@��@|�@;d@+@�@~��@~ȴ@~ȴ@~��@~v�@~V@~$�@}�@}�T@}��@}��@}`B@|��@|��@|I�@{��@{33@{@z�H@z��@z�\@zn�@zM�@zM�@zM�@z�@y��@y�7@yG�@y%@x�`@x�`@x�`@x�`@x�`@x��@x�@xQ�@xQ�@x �@w�@w�@w�@w�;@w��@w��@w�w@wl�@wK�@v��@v�@v�R@v��@v�+@vv�@v@u�-@u�h@u�h@u�h@u�@up�@uO�@u/@t��@t�/@tz�@s��@sC�@so@s@r�H@r��@r��@r�!@r�!@r�!@r�!@r�!@r�!@r�\@r-@q��@q��@q�^@q�^@q�^@q�7@qhs@qG�@qX@qX@qhs@qX@q%@pĜ@p�u@pbN@pQ�@p1'@pb@o�@o�@o�@o�@o\)@o\)@o\)@o\)@o\)@o\)@oK�@o�@n��@nȴ@nȴ@n�R@n��@n��@n�+@n�+@nV@nV@n$�@m�T@m�@mp�@mO�@m/@m�@l�@l�/@l��@l�j@l�j@l�@l��@l�D@l�D@lj@l1@k��@kƨ@k�F@k�F@k�F@k��@k�@kt�@kdZ@kC�@k"�@ko@j�@j��@j�!@j��@j��@j��@j��@j�\@j~�@jn�@jn�@jn�@j^5@j^5@j^5@jM�@j=q@j-@i��@i��@i�@i�#@i��@i��@i��@ix�@ihs@ihs@iX@iX@iX@iG�@i7L@i7L@i&�@h��@h�`@h�`@h�`@h�`@h�`@h�`@h��@h��@hĜ@hĜ@hĜ@h��@hĜ@h�9@h�9@h�9@h��@h�u@h�u@hr�@hbN@hbN@hbN@hbN@hbN@hr�@hr�@hr�@hr�@h�@h�@h�@h�@h�@h�`@h��@h�`@h�`@hĜ@h�9@h�u@hr�@hbN@hQ�@hQ�@hQ�@hQ�@hA�@hA�@hA�@h1'@h1'@h �@hb@hb@h �@h �@h �@h �@hb@hb@h �@hb@hb@hb@h  @hb@h �@h �@hb@h  @g�@g�;@g�;@g��@g�w@g�w@g�@g�@g��@g��@g�P@g�P@g|�@gl�@gK�@gK�@g;d@g;d@g+@g;d@g;d@g+@g�@g
=@f��@f��@f�y@f�y@fȴ@f�R@f��@f��@f�+@f�+@fv�@fff@fV@fE�@f5?@f$�@f{@f@e�T@e��@e@e�-@e�-@e��@e�h@e�@ep�@ep�@e`B@eO�@e?}@e?}@e?}@e/@e/@e/@e�@eV@eV@eV@d��G�O�@d�@d9X@c"�@a��@_�;@_+@_
=@^ff@_�P@`  @co@d(�@f��@h  @g\)@f@ct�@a��@`��@]�T@Z~�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBcTBcTBcTBcTBcTBcTBcTBcTBcTBcTBcTBcTBcTBcTBcTBbNBcTBcTBbNBcTBbNBbNBbNBbNBbNBbNBbNBbNBbNBcTBcTBcTBcTBcTBbNBcTBdZBgmBe`BiyBiyBiyBiyBiyBiyBiyBiyBiyBiyBiyBjBjBjBjBjBjBjBjBjBjBjBjBjBjBjBjBjBjBjBjBiyBiyBiyBiyBiyBiyBiyBiyBiyBjBjBjBjBjBjBk�Bl�Bm�Bn�Bp�Br�Bt�Bx�B~�B�B�+B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�{B�{B�{B�{B��B��B��B��B�{B�{B�{B�{B�{B�{B�{B�{B�{B�uB�uB�uB�uB�uB�uB�uB�uB�uB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�bB�hB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�\B�\B�\B�\B�\B�\B�bB�bB�bB�bB�bB�bB�bB�bB�bB�hB�hB�oB�oB�oB�hB�hB�hB�hB�bB�hB�hB�hB�hB�hB�hB�hB�hB�bB�bB�bB�bB�bB�bB�hB�hB�bB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�PB�VB�VB�PB�PB�PB�DB�+B�+B�+B�%B�7B�PB��B��B��B�!B�FB�XB�^B�XB�^B�jB�d1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 Bb�Bb�Bb�Bb�Bb�Bb�Bb�Bb�Bb�Bb�Bb�Bb�Bb�Bb�Bb�Ba�Bb�Bb�Ba�Bb�Ba�Ba�Ba�Ba�Ba�Ba�Ba�Ba�Ba�Bb�Bb�Bb�Bb�Bb�Ba�Bb�Bc�Bf�Bd�Bh�Bh�Bh�Bh�Bh�Bh�Bh�Bh�Bh�Bh�Bh�BjBjBjBjBjBjBjBjBj BjBj BjBjBjBjBjBjBjBjBjBh�Bh�Bh�Bh�Bh�Bh�Bh�Bh�Bh�BjBjBjBjBjBjBkBlBmBnBp(Br4Bt?Bx[B~|B��B��B�B�<B�6B�AB�HB�LB�BB�BB�CB�DB�CB�=B�<B�;B�;B�;B�<B�<B�<B�<B�6B�9B�=B�=B�:B�:B�:B�7B�9B�9B�6B�6B�6B�6B�6B�6B�6B�6B�7B�6B�8B�6B�8B�7B�8B�2B�2B�.B�.B�1B�/B�-B�(B�+B�,B�)B�*B�)B�'B�*B�+B�0B�0B�*B�*B�1B�1B�/B�1B�1B�/B�1B�,B�,B�,B�,B�,B�.B�#B�'B�,B�.B�,B�$B�)B�$B�*B�,B�$B�,B�,B�,B�*B�,B�%B�)B�+B�#B�#B�*B�+B�+B�*B�-B�*B�*B�+B�0B�/B�0B�/B�/B�0B�0B�/B�0B�/B�/B�+B�,B�*B�(B�*B�.B�*B�,B�*B�,B�+B�,B�+B�*B�%B�%B�%B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B�B��B�B�B�B�B�B�B�B��B��B��B�B� B�B� B� B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��G�O�B��B��B��B��B��B��B��B��B��B��B�B�8B��B��B��B��B��B��B��B��B��1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt( sw_cndr(PSAL,TEMP,PRES), TEMP, PRES_ADJUSTED )                                                                                                                                                                                         dP =-0.88 dbar.                                                                                                                                                                                                                                                 none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressures adjusted by using pressure offset at the sea surface. The quoted error is manufacturer specified accuracy in dbar.                                                                                                                                    The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   No significant salinity drift detected. The quoted error is max[0.01, 1xOWC uncertainty] in PSS-78.                                                                                                                                                             201904051344002019040513440020190405134400  AO  ARCAADJP                                                                    20181006004027    IP                G�O�G�O�G�O�                AO  ARGQQCPL                                                                    20181006004027  QCP$                G�O�G�O�G�O�5F03E           AO  ARGQQCPL                                                                    20181006004027  QCF$                G�O�G�O�G�O�8000            UW  ARSQUWQC    WOD & nearby Argo as visual check                               20190405134400  IP                  G�O�G�O�G�O�                