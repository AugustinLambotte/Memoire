CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   N_CALIB       	N_HISTORY             
   title         Argo float vertical profile    institution       AOML   source        
Argo float     history       2018-10-06T00:40:29Z creation      
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
_FillValue                    �PArgo profile    3.1 1.2 19500101000000  20181006004029  20190405134402  5904772 US ARGO PROJECT                                                 STEPHEN RISER                                                   PRES            TEMP            PSAL               CA   AO  6627                            2C  D   APEX                            7392                            062915                          846 @�R8�} 1   @�R9W:�@N�5?|��@O\(�1   GPS     Primary sampling: mixed [deeper than nominal 985dbar: discrete; nominal 985dbar to surface: 2dbar-bin averaged]                                                                                                                                                    CA   B   B   @���@�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(�C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D#��D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:�fD;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� D`��Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg�fDh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Dt��Dy� D�fD�FfD�s3D���D�3D�I�D�y�D���D�  D�9�D�l�D���D��D�I�D�P D�� D�  D�FfD�s3D��f111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @�  @�ffA33A+33AK33Ak33A���A���A���A���Ař�Aՙ�A噚A���B��B
��B��B��B"��B*��B2��B:��BB��BJ��BR��BZ��Bb��Bj��Br��Bz��B�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffB�ffC �3C�3C�3C�3C�3C
�3C�3C�3C�3C�3C�3C�3C�3C�3C�3C�3C �3C"�3C$�3C&�3C(��C*�3C,�3C.�3C0�3C2�3C4�3C6�3C8�3C:�3C<�3C>�3C@�3CB�3CD�3CF�3CH�3CJ�3CL�3CN�3CP�3CR�3CT�3CV�3CX�3CZ�3C\�3C^�3C`�3Cb�3Cd�3Cf�3Ch�3Cj�3Cl�3Cn�3Cp�3Cr�3Ct�3Cv�3Cx�3Cz�3C|�3C~�3C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�C�Y�D ,�D ��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D	,�D	��D
,�D
��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D,�D��D ,�D ��D!,�D!��D",�D"��D#,�D#��D$&gD$��D%,�D%��D&,�D&��D',�D'��D(,�D(��D),�D)��D*,�D*��D+,�D+��D,,�D,��D-,�D-��D.,�D.��D/,�D/��D0,�D0��D1,�D1��D2,�D2��D3,�D3��D4,�D4��D5,�D5��D6,�D6��D7,�D7��D8,�D8��D9,�D9��D:,�D:�3D;,�D;��D<,�D<��D=,�D=��D>,�D>��D?,�D?��D@,�D@��DA,�DA��DB,�DB��DC,�DC��DD,�DD��DE,�DE��DF,�DF��DG,�DG��DH,�DH��DI,�DI��DJ,�DJ��DK,�DK��DL,�DL��DM,�DM��DN,�DN��DO,�DO��DP,�DP��DQ,�DQ��DR,�DR��DS,�DS��DT,�DT��DU,�DU��DV,�DV��DW,�DW��DX,�DX��DY,�DY��DZ,�DZ��D[,�D[��D\,�D\��D],�D]��D^,�D^��D_,�D_��D`,�D`��Da&gDa��Db,�Db��Dc,�Dc��Dd,�Dd��De,�De��Df,�Df��Dg,�Dg�3Dh,�Dh��Di,�Di��Dj,�Dj��Dk,�Dk��Dl,�Dl��Dm,�Dm��Dn,�Dn��Do,�Do��Dp,�Dp��Dq,�Dq��Dr,�Dr��Ds,�Ds��Dt,�Dt��Dt��Dy��D��D�\�D���D��3D��D�` D�� D��3D�fD�P D��3D��3D�#3D�` D�ffD��fD�fD�\�D�D���111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���@���@��^@��^@��-@��7@�x�@��-@��-@���@���@�@�@���@�@�@�@�@���@��#@��@��@��@��@���@��@���@��-@���@��7@��h@��@��7@��h@��7@��@��@��@��@��@��@��@��@��@�p�@�x�@�x�@��@��7@��@��@��@��@��@�G�@�?}@�&�@���@�Ĝ@���@�z�@�z�@�Z@�1'@���@���@���@��F@�;d@���@���@��!@���@��@���@���@��y@�ȴ@���@���@��\@�$�@�O�@�Ĝ@��D@�Q�@��9@�1'@~v�@}/@{�@{o@z�\@zM�@z=q@z-@y��@y��@zJ@zJ@y�#@y��@y�^@y��@y�7@yG�@yhs@y%@x�9@x��@x�9@y��@y��@y�@z�H@|j@}/@}�-@}/@|��@{�@{t�@z�H@y�#@yhs@x��@x�`@x�`@x�`@x  @w�;@w
=@v��@v�R@vV@v$�@v@v@u�T@u��@u��@u�@v@v{@v5?@v5?@vE�@vV@vff@vv�@vff@v��@v�R@v�+@vv�@vV@vV@vff@v@v@v{@v$�@vff@v��@v��@v��@v�R@vȴ@vȴ@vȴ@vȴ@v��@vȴ@v�y@w+@w;d@w\)@wK�@wK�@w;d@wK�@wl�@w\)@v�y@v�R@v��@w�@wK�@wl�@w|�@w�P@w��@w�@w�@x��@z�@z��@z��@{dZ@{ƨ@{�m@{�
@{��@{t�@{t�@{�
@|9X@|9X@|9X@|�@{�m@{��@{C�@{C�@{"�@zn�@zJ@y�@y��@zJ@yX@x��@xbN@xA�@xr�@xr�@x��@x�`@x�`@x��@x�u@xQ�@xA�@x  @w�w@w�P@w|�@wl�@w\)@wK�@w\)@w\)@w|�@w|�@w;d@w+@w+@v�R@v5?@u�h@u�@u?}@u@v{@u�h@t�@s�m@s�m@s�
@s��@s��@st�@s�@s�@s��@s��@r^5@qX@qX@q��@q��@q��@q��@q��@r=q@r=q@r=q@rn�@r~�@r~�@rM�@q��@q�@q��@q��@q�#@q�#@q��@q&�@pĜ@p��@o�@o;d@n�R@nff@nE�@m��@m�@m��@m��@mp�@m�@l�@l�j@l�@l9X@l�@l1@k��@kƨ@k��@kt�@kt�@kdZ@kC�@j�!@jM�@i�#@ihs@iG�@iG�@i&�@i&�@i7L@i7L@iX@ihs@h�`@h�u@h�@h�u@h�u@h��@h�u@h�u@h�@hr�@h �@g��@g�;@g��@g�;@g�@g�@g��@g�w@g��@g��@g�P@gl�@g\)@g+@g;d@g;d@g+@g+@g+@g+@g�@g
=@g
=@g
=@g
=@g
=@f��@f�y@f�@fȴ@f��@f��@f�+@fv�@fv�@fv�@fv�@f5?@e�@e�-@e�-@e��@e�h@e�h@e�@e�@e�@ep�@e?}@e�@e�@eV@d��@d��@d�@d�/@d�D@dz�@dj@d(�@c��@c�
@c�
@cƨ@cƨ@cƨ@c��@c�@ct�@ct�@ct�@cdZ@cS�@cC�@c33@c"�@c"�@c"�@c"�@c"�@co@c@b�@b��@bn�@bn�@b^5@b-@bJ@a��@a�@a�^@a��@a��@a�^@a�7@ax�@aX@a7L@a7L@a&�@a%@`��@`��@`�`@`��@`�`@a%@a%@`��@`Ĝ@`Ĝ@`��@`��@`�`@`��@`��@`��@`Ĝ@`Ĝ@`Ĝ@`�9@`�9@`�u@`bN@`bN@`A�@`  @_�;@_�P@_�P@_|�@_l�@_l�@_K�@_�@_�@_�@_�@^��@^�y@^�@^�R@^��@^��@^��@^��@^�+@^�+@^�+@^ff@^$�@]`B@\(�@[�@Z�!@Z��@Z�@^�y@a7L@d��@g+@i%@i%@f��@d�j@a��@^�R@\�@Y��@XbN@Up�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  @���@���@��^@��^@��-@��7@�x�@��-@��-@���@���@�@�@���@�@�@�@�@���@��#@��@��@��@��@���@��@���@��-@���@��7@��h@��@��7@��h@��7@��@��@��@��@��@��@��@��@��@�p�@�x�@�x�@��@��7@��@��@��@��@��@�G�@�?}@�&�@���@�Ĝ@���@�z�@�z�@�Z@�1'@���@���@���@��F@�;d@���@���@��!@���@��@���@���@��y@�ȴ@���@���@��\@�$�@�O�@�Ĝ@��D@�Q�@��9@�1'@~v�@}/@{�@{o@z�\@zM�@z=q@z-@y��@y��@zJ@zJ@y�#@y��@y�^@y��@y�7@yG�@yhs@y%@x�9@x��@x�9@y��@y��@y�@z�H@|j@}/@}�-@}/@|��@{�@{t�@z�H@y�#@yhs@x��@x�`@x�`@x�`@x  @w�;@w
=@v��@v�R@vV@v$�@v@v@u�T@u��@u��@u�@v@v{@v5?@v5?@vE�@vV@vff@vv�@vff@v��@v�R@v�+@vv�@vV@vV@vff@v@v@v{@v$�@vff@v��@v��@v��@v�R@vȴ@vȴ@vȴ@vȴ@v��@vȴ@v�y@w+@w;d@w\)@wK�@wK�@w;d@wK�@wl�@w\)@v�y@v�R@v��@w�@wK�@wl�@w|�@w�P@w��@w�@w�@x��@z�@z��@z��@{dZ@{ƨ@{�m@{�
@{��@{t�@{t�@{�
@|9X@|9X@|9X@|�@{�m@{��@{C�@{C�@{"�@zn�@zJ@y�@y��@zJ@yX@x��@xbN@xA�@xr�@xr�@x��@x�`@x�`@x��@x�u@xQ�@xA�@x  @w�w@w�P@w|�@wl�@w\)@wK�@w\)@w\)@w|�@w|�@w;d@w+@w+@v�R@v5?@u�h@u�@u?}@u@v{@u�h@t�@s�m@s�m@s�
@s��@s��@st�@s�@s�@s��@s��@r^5@qX@qX@q��@q��@q��@q��@q��@r=q@r=q@r=q@rn�@r~�@r~�@rM�@q��@q�@q��@q��@q�#@q�#@q��@q&�@pĜ@p��@o�@o;d@n�R@nff@nE�@m��@m�@m��@m��@mp�@m�@l�@l�j@l�@l9X@l�@l1@k��@kƨ@k��@kt�@kt�@kdZ@kC�@j�!@jM�@i�#@ihs@iG�@iG�@i&�@i&�@i7L@i7L@iX@ihs@h�`@h�u@h�@h�u@h�u@h��@h�u@h�u@h�@hr�@h �@g��@g�;@g��@g�;@g�@g�@g��@g�w@g��@g��@g�P@gl�@g\)@g+@g;d@g;d@g+@g+@g+@g+@g�@g
=@g
=@g
=@g
=@g
=@f��@f�y@f�@fȴ@f��@f��@f�+@fv�@fv�@fv�@fv�@f5?@e�@e�-@e�-@e��@e�h@e�h@e�@e�@e�@ep�@e?}@e�@e�@eV@d��@d��@d�@d�/@d�D@dz�@dj@d(�@c��@c�
@c�
@cƨ@cƨ@cƨ@c��@c�@ct�@ct�@ct�@cdZ@cS�@cC�@c33@c"�@c"�@c"�@c"�@c"�@co@c@b�@b��@bn�@bn�@b^5@b-@bJ@a��@a�@a�^@a��@a��@a�^@a�7@ax�@aX@a7L@a7L@a&�@a%@`��@`��@`�`@`��@`�`@a%@a%@`��@`Ĝ@`Ĝ@`��@`��@`�`@`��@`��@`��@`Ĝ@`Ĝ@`Ĝ@`�9@`�9@`�u@`bN@`bN@`A�@`  @_�;@_�P@_�P@_|�@_l�@_l�@_K�@_�@_�@_�@_�@^��@^�y@^�@^�R@^��@^��@^��@^��@^�+@^�+@^�+G�O�@^$�@]`B@\(�@[�@Z�!@Z��@Z�@^�y@a7L@d��@g+@i%@i%@f��@d�j@a��@^�R@\�@Y��@XbN@Up�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�oB�oB�oB�oB�oB�uB�uB�oB�oB�oB�oB�hB�hB�hB�hB�hB�hB�hB�hB�oB�oB�oB�oB�oB�oB�hB�hB�hB�hB�hB�bB�hB�bB�bB�bB�bB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�oB�oB�oB�oB�uB�uB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�uB�hB�\B�VB�PB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�DB�DB�DB�JB�DB�DB�JB�JB�VB�\B�\B�hB�{B��B��B��B��B�{B�uB�oB�bB�bB�\B�\B�\B�\B�VB�VB�JB�JB�DB�DB�DB�=B�=B�=B�=B�=B�DB�DB�DB�JB�JB�JB�JB�JB�PB�PB�PB�VB�VB�PB�PB�PB�PB�PB�PB�PB�PB�VB�\B�\B�\B�\B�\B�bB�bB�bB�bB�bB�hB�hB�hB�oB�oB�oB�oB�oB�oB�oB�oB�oB�uB�uB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B��B��B��B��B��B�{B�{B�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�oB�oB�oB�oB�oB�oB�oB�oB�oB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�\B�\B�\B�\B�\B�\B�\B�\B�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�PB�PB�PB�PB�PB�PB�PB�PB�JB�JB�JB�DB�DB�DB�DB�DB�DB�DB�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�7B�7B�7B�7B�1B�1B�1B�1B�1B�1B�+B�+B�1B�1B�+B�+B�+B�+B�+B�%B�%B�%B�%B�%B�%B�%B�+B�+B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B~�B}�B|�B|�B� B�DB�{B��B�B�?B�XB�jB�qB�qB�jB�jB�jB�jB�q111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111114  B�
B�B�
B�B�B�B�B�B�B�
B�B�B�B�B�B�B�B�B�B�B�
B�
B�
B�
B�B�B�B�B�B�B� B�B��B��B� B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�	B�B�B�B�B�	B�	B�
B�B�B�
B�B�>B�BB�BB�IB�GB�IB�NB�OB�OB�PB�OB�bB�bB�bB�ZB�ZB�TB�DB�6B�6B�0B�6B�"B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�"B�+B�,B�%B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�
B�
B�
B�B�
B�
B�B�
B�B�B�B�B�B�B�&B�)B�(B�=B�OB�ZB�YB�lB�pB�uB�uB�{B�zB�}B��B��B��B��B��B��B��B��B��B�B�zB�xB�zB��B��B�xB�yB�lB�nB�vB�vB�yB�yB��B��B��B��B�{B�{B�|B�{B�{B�{B�uB�|B�yB��B��B��B��B��B��B��B��B�|B�wB�|B��B��B��B�|B�vB�vB�vB�nB�nB�pB�vB�nB�vB�pB�gB�]B�]B�cB�hB�hB�hB�hB�qB�qB�oB�qB�oB�oB�qB�oB�qB�oB�oB�qB�jB�jB�eB�cB�\B�VB�PB�KB�FB�FB�?B�FB�?B�?B�?B�8B�4B�1B�2B�+B�+B�+B�+B�'B�B�'B�)B�'B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�
B�
B�B�
B�
B�B�B�B�B�B�B�B�B�B�B� B�B�B�B�B�B� B�B� B�B�B�B� B� B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��G�O�B��B��B~�B}�B|�B|�B�B��B�#B�rB��B��B��B�B�B�B�B�B�B�G�O�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111114  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
G�O�PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt( sw_cndr(PSAL,TEMP,PRES), TEMP, PRES_ADJUSTED )                                                                                                                                                                                         dP =-0.7 dbar.                                                                                                                                                                                                                                                  none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressures adjusted by using pressure offset at the sea surface. The quoted error is manufacturer specified accuracy in dbar.                                                                                                                                    The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   No significant salinity drift detected. The quoted error is max[0.01, 1xOWC uncertainty] in PSS-78.                                                                                                                                                             201904051344022019040513440220190405134402  AO  ARCAADJP                                                                    20181006004029    IP                G�O�G�O�G�O�                AO  ARGQQCPL                                                                    20181006004029  QCP$                G�O�G�O�G�O�5F03E           AO  ARGQQCPL                                                                    20181006004029  QCF$                G�O�G�O�G�O�8000            UW  ARSQUWQC    WOD & nearby Argo as visual check                               20190405134402  IP                  G�O�G�O�G�O�                