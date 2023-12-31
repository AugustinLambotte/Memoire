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
_FillValue                    �PArgo profile    3.1 1.2 19500101000000  20181006004029  20190405134402  5904772 US ARGO PROJECT                                                 STEPHEN RISER                                                   PRES            TEMP            PSAL               FA   AO  6627                            2C  D   APEX                            7392                            062915                          846 @�Y��۲1   @�Y�<���@O��7Kƨ�@��+1   GPS     Primary sampling: mixed [deeper than nominal 985dbar: discrete; nominal 985dbar to surface: 2dbar-bin averaged]                                                                                                                                                    FA   B   B   @�33@�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Dt� Dy��D�	�D�L�D��3D���D� D�@ D���D�ɚD�3D�I�D�p D��3D� D�@ D�c3D���D�fD�9�D� D��f111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@�Q�A(�A,(�AL(�Al(�A�{A�{A�{A�{A�{A�{A�{A�{B
=B
=B
=B
=B#
=B+
=B3
=B;
=BC
=BK
=BS
=B[
=Bc
=Bk
=Bs
=B{
=B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BŅBɅBͅBхBՅBمB݅B�B�B�B�B�B��B��B��C CCCCC
CCCCCCCCCCC C"C$C&C(C*C,C.C0C2C4C6C8C:C<C>C@CBCDCFCHCJCLCNCPCRCTCVCXCZC\C^C`CbCdCfChCjClCnCpCrCtCvCxCzC|C~C�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�nC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�T{C�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHD 0�D ��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D	0�D	��D
0�D
��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D 0�D ��D!0�D!��D"0�D"��D#0�D#��D$0�D$��D%0�D%��D&0�D&��D'0�D'��D(0�D(��D)0�D)��D*0�D*��D+0�D+��D,0�D,��D-0�D-��D.0�D.��D/0�D/��D00�D0��D10�D1��D20�D2��D30�D3��D40�D4��D50�D5��D60�D6��D70�D7��D80�D8��D90�D9��D:0�D:��D;0�D;��D<0�D<��D=0�D=��D>0�D>��D?0�D?��D@0�D@��DA0�DA��DB0�DB��DC0�DC��DD0�DD��DE0�DE��DF0�DF��DG0�DG��DH0�DH��DI0�DI��DJ0�DJ��DK0�DK��DL0�DL��DM0�DM��DN0�DN��DO0�DO��DP0�DP��DQ0�DQ��DR0�DR��DS0�DS��DT0�DT��DU0�DU��DV0�DV��DW0�DW��DX0�DX��DY0�DY��DZ0�DZ��D[0�D[��D\0�D\��D]0�D]��D^0�D^��D_0�D_��D`0�D`��Da0�Da��Db0�Db��Dc0�Dc��Dd0�Dd��De0�De��Df0�Df��Dg0�Dg��Dh0�Dh��Di0�Di��Dj0�Dj��Dk0�Dk��Dl0�Dl��Dm0�Dm��Dn0�Dn��Do0�Do��Dp0�Dp��Dq0�Dq��Dr0�Dr��Ds0�Ds��Dt0�Dt��Du�Dy�>D�!�D�eD���D��D�(RD�XRD��D���D�+�D�a�D��RD�ۅD�(RD�XRD�{�D��D��D�Q�D�RD���111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�Z@�I�@�A�@�1'@�1'@�b@�  @�  @�1@�  @���@��@��@��m@��m@��
@��
@��w@��F@��F@���@���@��@�C�@�33@�o@��@�~�@��^@��/@� �@���@�7L@��w@�v�@�^5@�E�@�5?@�$�@�$�@�$�@�$�@��@��@��@��T@��^@�@��^@���@���@��@�G�@�G�@�G�@�?}@�/@�&�@��@�%@��/@���@��@�j@�(�@�ƨ@���@�\)@�+@��H@��R@�n�@�-@��@���@���@��h@�hs@�?}@��@��/@��j@��j@��j@��@��D@��@�z�@�z�@�r�@�Z@�I�@�A�@�A�@�9X@�1'@� �@�b@�1@�  @��@���@���@�t�@�dZ@�K�@�@���@��R@��R@��R@���@���@�ȴ@���@��y@��H@��\@�"�@��u@��u@�dZ@�@��`@�j@�Q�@�9X@���@��w@��P@�;d@��!@�v�@��@��#@���@�7L@�%@��@�bN@�b@��@��F@���@��@���@��@�33@��@�33@��y@��@���@��@���@�hs@�?}@�O�@�?}@�x�@�{@�^5@���@��!@��\@��\@���@���@�@�o@��@�r�@�&�@��@��F@��
@�b@� �@���@�dZ@���@�M�@��@��-@�G�@��j@� �@�ƨ@�t�@�K�@�"�@��H@��@���@�v�@�M�@�@���@���@�`B@�hs@�E�@��m@�(�@���@�dZ@���@��R@���@�n�@�E�@�-@��-@�x�@��@�p�@�O�@�%@�bN@���@��@�K�@�o@�ȴ@��@��7@�p�@�hs@�O�@�?}@�7L@���@��D@�bN@�9X@��@�1@���@���@�t�@�l�@�S�@�33@��H@���@���@��+@�n�@�V@���@���@���@��/@�/@���@�I�@�I�@�A�@�(�@� �@��@�1@���@��m@���@�|�@�S�@�;d@�o@��y@�n�@�=q@�J@�$�@�{@��#@�@��^@��-@��^@��T@���@�-@�J@���@�O�@��@���@��j@��9@���@��u@�Z@��@��;@���@�K�@��@�+@���@���@���@��!@��!@��!@���@���@��+@�^5@�{@���@��-@��h@�x�@�p�@�`B@�r�@�dZ@��R@�M�@���@���@�Q�@�(�@��@~�y@~�@K�@~��@~�+@}�T@|��@|�/@|��@|9X@{�
@{33@z�@{@z�@{o@{33@{dZ@{��@{dZ@z��@z^5@z=q@y��@y�#@y��@y��@yx�@yX@yX@yG�@y7L@y&�@y�@x��@x��@x�@xbN@xQ�@xQ�@x �@x  @w�@w�@w�@w�;@w��@w��@w�w@w��@w�w@w��@w�P@w�P@w|�@wl�@wK�@w;d@w+@w+@w
=@v�y@v�@v�@vȴ@v�R@v�R@v�R@v�R@v��@v�+@v�y@w�P@w�@xr�@y�^@y��@y��@yhs@yhs@yX@yx�@y��@zJ@z^5@{o@{o@z~�@x��@xA�@x�9@xĜ@xĜ@x�@vE�@t��@sƨ@sƨ@s�F@s��@tZ@vE�@vv�@w�@w\)@w|�@w\)@wK�@wK�@w�@w�@x  @w�;@w�@wK�@w
=@v��@v��@u��@u?}@uV@t��@r~�@q��@qG�@p�`@pĜ@pĜ@pĜ@pĜ@p��@pĜ@q�@qG�@p��@p��@p��@p�9@pĜ@p�9@pĜ@p��@p��@pĜ@p��@p��@pĜ@p��@q�@qX@q�7@q�#@rJ@rJ@r�@r�@r-@rJ@q�@q�^@q��@qx�@q7L@q&�@o�@lI�@i�^@g��@f�+@c�m@cƨ@b-@`��@^ff@\��@\1@[dZ@Z=q@a��@bn�@h�@kƨ@l(�@i&�@f�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  @�Z@�I�@�A�@�1'@�1'@�b@�  @�  @�1@�  @���@��@��@��m@��m@��
@��
@��w@��F@��F@���@���@��@�C�@�33@�o@��@�~�@��^@��/@� �@���@�7L@��w@�v�@�^5@�E�@�5?@�$�@�$�@�$�@�$�@��@��@��@��T@��^@�@��^@���@���@��@�G�@�G�@�G�@�?}@�/@�&�@��@�%@��/@���@��@�j@�(�@�ƨ@���@�\)@�+@��H@��R@�n�@�-@��@���@���@��h@�hs@�?}@��@��/@��j@��j@��j@��@��D@��@�z�@�z�@�r�@�Z@�I�@�A�@�A�@�9X@�1'@� �@�b@�1@�  @��@���@���@�t�@�dZ@�K�@�@���@��R@��R@��R@���@���@�ȴ@���@��y@��H@��\@�"�@��u@��u@�dZ@�@��`@�j@�Q�@�9X@���@��w@��P@�;d@��!@�v�@��@��#@���@�7L@�%@��@�bN@�b@��@��F@���@��@���@��@�33@��@�33@��y@��@���@��@���@�hs@�?}@�O�@�?}@�x�@�{@�^5@���@��!@��\@��\@���@���@�@�o@��@�r�@�&�@��@��F@��
@�b@� �@���@�dZ@���@�M�@��@��-@�G�@��j@� �@�ƨ@�t�@�K�@�"�@��H@��@���@�v�@�M�@�@���@���@�`B@�hs@�E�@��m@�(�@���@�dZ@���@��R@���@�n�@�E�@�-@��-@�x�@��@�p�@�O�@�%@�bN@���@��@�K�@�o@�ȴ@��@��7@�p�@�hs@�O�@�?}@�7L@���@��D@�bN@�9X@��@�1@���@���@�t�@�l�@�S�@�33@��H@���@���@��+@�n�@�V@���@���@���@��/@�/@���@�I�@�I�@�A�@�(�@� �@��@�1@���@��m@���@�|�@�S�@�;d@�o@��y@�n�@�=q@�J@�$�@�{@��#@�@��^@��-@��^@��T@���@�-@�J@���@�O�@��@���@��j@��9@���@��u@�Z@��@��;@���@�K�@��@�+@���@���@���@��!@��!@��!@���@���@��+@�^5@�{@���@��-@��h@�x�@�p�@�`B@�r�@�dZ@��R@�M�@���@���@�Q�@�(�@��@~�y@~�@K�@~��@~�+@}�T@|��@|�/@|��@|9X@{�
@{33@z�@{@z�@{o@{33@{dZ@{��@{dZ@z��@z^5@z=q@y��@y�#@y��@y��@yx�@yX@yX@yG�@y7L@y&�@y�@x��@x��@x�@xbN@xQ�@xQ�@x �@x  @w�@w�@w�@w�;@w��@w��@w�w@w��@w�w@w��@w�P@w�P@w|�@wl�@wK�@w;d@w+@w+@w
=@v�y@v�@v�@vȴ@v�R@v�R@v�R@v�R@v��@v�+@v�y@w�P@w�@xr�@y�^@y��@y��@yhs@yhs@yX@yx�@y��@zJ@z^5@{o@{o@z~�@x��@xA�@x�9@xĜ@xĜ@x�@vE�@t��@sƨ@sƨ@s�F@s��@tZ@vE�@vv�@w�@w\)@w|�@w\)@wK�@wK�@w�@w�@x  @w�;@w�@wK�@w
=@v��@v��@u��@u?}@uV@t��@r~�@q��@qG�@p�`@pĜ@pĜ@pĜ@pĜ@p��@pĜ@q�@qG�@p��@p��@p��@p�9@pĜ@p�9@pĜ@p��@p��@pĜ@p��@p��@pĜ@p��@q�@qX@q�7@q�#@rJ@rJ@r�@r�@r-@rJ@q�@q�^@q��@qx�@q7LG�O�@o�@lI�@i�^@g��@f�+@c�m@cƨ@b-@`��@^ff@\��@\1@[dZ@Z=q@a��@bn�@h�@kƨ@l(�@i&�@f�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BǮBŢBĜBÖBÖBÖBÖBÖBÖBÖBÖBÖBBB��B��B��B��B��B��B��B�}B��B��B�}B�}B�}B�}B�wB�wB�qB�jB�dB�^B�XB�XB�RB�RB�LB�FB�FB�?B�?B�?B�?B�?B�?B�9B�9B�9B�9B�9B�9B�9B�3B�3B�3B�3B�3B�3B�3B�3B�3B�3B�3B�3B�3B�3B�-B�-B�'B�'B�'B�!B�!B�!B�!B�!B�!B�!B�!B�!B�'B�-B�-B�-B�FB�qB�qB�^B�?B�'B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�'B�B�B�!B�3B�3B�3B�-B�!B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�!B�RB�^B�XB�RB�LB�FB�FB�?B�?B�9B�9B�9B�?B�?B�9B�3B�-B�!B�!B�!B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�!B�'B�!B�B�B�B�B�!B�!B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B��B��B�{B�uB�uB�uB�uB�uB�uB�{B�{B��B�{B�uB�uB�uB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�hB�hB�hB�hB�hB�hB�hB�hB�oB�oB�oB�oB�oB�oB�oB�oB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�{B�{B�{B�{B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�uB�uB�uB�oB�uB�{B�{B�{B��B��B��B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�oB�\B�PB�1B�7B�%B�B~�B|�B{�By�Bv�B�VB�oB��B�9B�qB�wB�}111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  BԖBԘBԖBԖBԘBԖBԖBԖBԖBԘBԖBԖBԖBԖBԖBԘBԘBԖBԖBԖBԖBԖBӏBӍBӏB҉B҇BцB�{B�nB�aB�XB�CB�9B�3B�.B�.B�-B�-B�-B�-B�+B�-B�.B�$B�$B�!B�B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B�B�zB�uB�pB�gB�cB�[B�UB�VB�UB�MB�UB�VB�VB�LB�MB�OB�KB�PB�JB�6B�2B�,B�+B�+B�+B�6B�BB�PB�VB�^B�bB�dB�uB�sB�uB�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�pB�mB�uB�yB�{B�vB�uB�wB�uB�tB�uB�uB�nB�pB�sB�vB�vB�sB�uB�sB�oB�qB�tB�|B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�iB�UB�QB�JB�2B�+B�*B�&B�B�B�%B�$B�B�B�B�B�B�B�
B�
B�
B�
B�B�B�B�B�B�B�B�B�B�B�B�	B�B�B�B�	B�B�B�B�B��B��B��B��B��B��B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�"B�%B�/B�9B�KB�QB�NB�IB�LB�KB�TB�RB�XB�]B�eB�eB�_B�NB�DB�OB�QB�QB�EB�1B�B�B�B�B�B�%B�EB�JB�RB�XB�VB�WB�WB�WB�`B�dB�dB�eB�`B�`B�YB�WB�OB�CB�CB�@B�2B�"B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B� B� B� B�B� B� B�B�'B�+B�1B�8B�BB�?B�FB�GB�FB�DB�FB�FB�FB�FB�@B�@G�O�B�<B�B�B��B��B��B��B��B��B~�B|�B{�ByxBvfB��B�B��B��B�B�B� 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt( sw_cndr(PSAL,TEMP,PRES), TEMP, PRES_ADJUSTED )                                                                                                                                                                                         dP =-0.76 dbar.                                                                                                                                                                                                                                                 none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressures adjusted by using pressure offset at the sea surface. The quoted error is manufacturer specified accuracy in dbar.                                                                                                                                    The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   No significant salinity drift detected. The quoted error is max[0.01, 1xOWC uncertainty] in PSS-78.                                                                                                                                                             201904051344022019040513440220190405134402  AO  ARCAADJP                                                                    20181006004029    IP                G�O�G�O�G�O�                AO  ARGQQCPL                                                                    20181006004029  QCP$                G�O�G�O�G�O�5F03E           AO  ARGQQCPL                                                                    20181006004029  QCF$                G�O�G�O�G�O�8000            UW  ARSQUWQC    WOD & nearby Argo as visual check                               20190405134402  IP                  G�O�G�O�G�O�                