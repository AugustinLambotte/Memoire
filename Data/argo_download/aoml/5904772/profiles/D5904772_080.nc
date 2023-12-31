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
_FillValue                    �PArgo profile    3.1 1.2 19500101000000  20181006004031  20190405134404  5904772 US ARGO PROJECT                                                 STEPHEN RISER                                                   PRES            TEMP            PSAL               PA   AO  6627                            2C  D   APEX                            7392                            062915                          846 @�r�
n1   @�r����@Oa���l��C;C��%1   GPS     Primary sampling: mixed [deeper than nominal 985dbar: discrete; nominal 985dbar to surface: 2dbar-bin averaged]                                                                                                                                                    PA   B   B   @�  @�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5fD5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Dt�fDy��D�fD�9�D�` D���D� D�,�D�� D��3D���D�,�D��3D��fD� D�6fD�vfD��fD��fD�I�D�c3D��f111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @�Q�@�Q�A(�A,(�AL(�Al(�A�{A�{A�{A�{A�{A�{A�{A�{B
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
CCCCCCCCCCC C"C$C&C(C*C,C.C0C2C4C6C8C:C<C>C@CBCDCFCHCJCLCNCPCRCTCVCXCZC\C^C`CbCdCfChCjClCnCpCrCtCvCxCzC|C~C�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHC�aHD 0�D ��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D	0�D	��D
0�D
��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D0�D��D 0�D ��D!0�D!��D"0�D"��D#0�D#��D$0�D$��D%0�D%��D&0�D&��D'0�D'��D(0�D(��D)0�D)��D*0�D*��D+0�D+��D,0�D,��D-0�D-��D.0�D.��D/0�D/��D00�D0��D10�D1��D20�D2��D30�D3��D40�D4��D57
D5��D60�D6��D70�D7��D80�D8��D90�D9��D:0�D:��D;0�D;��D<0�D<��D=0�D=��D>0�D>��D?0�D?��D@0�D@��DA0�DA��DB0�DB��DC0�DC��DD0�DD��DE0�DE��DF0�DF��DG0�DG��DH0�DH��DI0�DI��DJ0�DJ��DK0�DK��DL0�DL��DM0�DM��DN0�DN��DO0�DO��DP0�DP��DQ0�DQ��DR0�DR��DS0�DS��DT0�DT��DU0�DU��DV0�DV��DW0�DW��DX0�DX��DY0�DY��DZ0�DZ��D[0�D[��D\0�D\��D]0�D]��D^0�D^��D_0�D_��D`0�D`��Da0�Da��Db0�Db��Dc0�Dc��Dd0�Dd��De0�De��Df0�Df��Dg0�Dg��Dh0�Dh��Di0�Di��Dj0�Dj��Dk0�Dk��Dl0�Dl��Dm0�Dm��Dn0�Dn��Do0�Do��Dp0�Dp��Dq0�Dq��Dr0�Dr��Ds0�Ds��Dt0�Dt��Du
Dy�qD�.�D�Q�D�xRD��D�(RD�ED��RD�ۅD�D�ED���D�޸D�(RD�N�Dڎ�D��D��D�a�D�{�D�޸111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��A1'A�A�FA��Al�A	�^AffA�A��@�I�@���@��@�-@�?}@�?}@��m@�Z@��@ȼj@�X@�V@�1'@��
@�dZ@��@��7@�V@��@�Z@� �@��m@�+@��\@��T@�`B@�A�@���@��H@���@��7@�x�@�x�@�X@��@���@�Q�@�  @��
@�1@�t�@���@�5?@��^@�&�@��@���@�dZ@��H@�=q@�@��h@�?}@���@�r�@��@���@���@�t�@�S�@�^5@��7@�O�@���@���@���@��@���@��@�;d@�33@��@��R@���@���@�^5@�/@�Ĝ@��u@�z�@��@�j@��@�^5@�-@��T@�&�@��@���@�@��@��@�@��7@�7L@���@�A�@��@�b@��;@��w@���@�t�@�;d@�
=@��H@���@�M�@�@���@���@���@�Z@�I�@���@�t�@��y@��\@�~�@�$�@��7@��@��h@�x�@�`B@��@��j@��@��@�  @��@�l�@�C�@���@�n�@�n�@�v�@��@��#@��@�-@�$�@�J@��@���@�X@���@�A�@��@�\)@�o@��+@�@���@�G�@�Ĝ@�I�@���@�;d@�"�@�@���@��!@��+@�5?@��@�x�@�X@�/@�V@���@���@��9@���@���@�j@�ƨ@��y@�^5@�=q@�-@�{@��@��#@�@���@��7@�hs@�Ĝ@��@��9@��j@���@���@��@�%@���@���@���@�Ĝ@��j@��@���@���@��u@�bN@�I�@�I�@��@�K�@�J@���@���@�%@��D@�(�@��@�|�@�dZ@�l�@��@��y@���@���@�V@�^5@�v�@�V@�{@�@�@�@��@��^@�hs@�&�@�%@��j@��D@�r�@�bN@�Z@�b@
=@~�y@~�@~�y@~�@~�y@+@\)@l�@��@�@�@�w@�@~@}��@}O�@|�j@{��@{C�@{C�@{33@{33@{o@z��@zM�@zJ@y��@yx�@yx�@y�7@y�7@y�7@yX@y7L@y7L@y7L@y7L@yG�@yG�@yG�@x�`@x�u@x�@x1'@xb@w��@w��@w�w@w�;@xbN@xr�@xbN@xbN@xQ�@xr�@x�u@x��@y7L@y%@x�@x1'@x  @w�@w�P@w+@v��@v�R@vff@vE�@v5?@v$�@v@v@u�@u�T@u��@u��@up�@t��@t�@t��@tz�@tj@tZ@t9X@t(�@t�@t1@sƨ@s��@st�@sS�@sS�@s"�@r�@r��@r��@r��@r��@r��@r�\@r�@q�7@q&�@p�9@p��@p��@p�9@p�u@pA�@o��@o�w@ol�@o;d@o
=@o
=@o
=@o
=@n��@n�y@n�y@o�@o�@o+@oK�@o\)@o\)@o�P@ol�@o+@p  @p  @o�;@o��@o|�@o+@o+@n�y@nȴ@n��@n�@n�@nȴ@n��@nV@m�T@m/@l�@l�/@m`B@m`B@mp�@m�-@m�-@m�@m/@m�@l�@l��@l�D@l�D@l��@m�@m�@m�@m�@l��@l�j@l�D@l�D@l��@l�j@l��@l�j@l�@l�@l�D@l�@k��@l1@l1@k��@k��@k��@k��@kt�@k@j^5@jn�@j=q@i�#@i�^@i��@i��@i�^@i�^@i��@iX@i��@i��@i��@iX@i&�@i&�@i%@h��@i7L@hĜ@hQ�@hQ�@h �@hb@h  @hb@hb@g�w@g|�@g|�@g|�@g\)@g;d@g;d@g;d@g+@f�y@f�R@f�R@f��@f�+@f�+@f�+@fV@f{@f@e�@e��@d�@d9X@aG�@`1'@^$�@_
=@`Q�@b�H@d�@f$�@f�R@f�+@dZ@`�u@\1@O�@J�H@J-@:�H@2�!@+t�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  A1'A�A�FA��Al�A	�^AffA�A��@�I�@���@��@�-@�?}@�?}@��m@�Z@��@ȼj@�X@�V@�1'@��
@�dZ@��@��7@�V@��@�Z@� �@��m@�+@��\@��T@�`B@�A�@���@��H@���@��7@�x�@�x�@�X@��@���@�Q�@�  @��
@�1@�t�@���@�5?@��^@�&�@��@���@�dZ@��H@�=q@�@��h@�?}@���@�r�@��@���@���@�t�@�S�@�^5@��7@�O�@���@���@���@��@���@��@�;d@�33@��@��R@���@���@�^5@�/@�Ĝ@��u@�z�@��@�j@��@�^5@�-@��T@�&�@��@���@�@��@��@�@��7@�7L@���@�A�@��@�b@��;@��w@���@�t�@�;d@�
=@��H@���@�M�@�@���@���@���@�Z@�I�@���@�t�@��y@��\@�~�@�$�@��7@��@��h@�x�@�`B@��@��j@��@��@�  @��@�l�@�C�@���@�n�@�n�@�v�@��@��#@��@�-@�$�@�J@��@���@�X@���@�A�@��@�\)@�o@��+@�@���@�G�@�Ĝ@�I�@���@�;d@�"�@�@���@��!@��+@�5?@��@�x�@�X@�/@�V@���@���@��9@���@���@�j@�ƨ@��y@�^5@�=q@�-@�{@��@��#@�@���@��7@�hs@�Ĝ@��@��9@��j@���@���@��@�%@���@���@���@�Ĝ@��j@��@���@���@��u@�bN@�I�@�I�@��@�K�@�J@���@���@�%@��D@�(�@��@�|�@�dZ@�l�@��@��y@���@���@�V@�^5@�v�@�V@�{@�@�@�@��@��^@�hs@�&�@�%@��j@��D@�r�@�bN@�Z@�b@
=@~�y@~�@~�y@~�@~�y@+@\)@l�@��@�@�@�w@�@~@}��@}O�@|�j@{��@{C�@{C�@{33@{33@{o@z��@zM�@zJ@y��@yx�@yx�@y�7@y�7@y�7@yX@y7L@y7L@y7L@y7L@yG�@yG�@yG�@x�`@x�u@x�@x1'@xb@w��@w��@w�w@w�;@xbN@xr�@xbN@xbN@xQ�@xr�@x�u@x��@y7L@y%@x�@x1'@x  @w�@w�P@w+@v��@v�R@vff@vE�@v5?@v$�@v@v@u�@u�T@u��@u��@up�@t��@t�@t��@tz�@tj@tZ@t9X@t(�@t�@t1@sƨ@s��@st�@sS�@sS�@s"�@r�@r��@r��@r��@r��@r��@r�\@r�@q�7@q&�@p�9@p��@p��@p�9@p�u@pA�@o��@o�w@ol�@o;d@o
=@o
=@o
=@o
=@n��@n�y@n�y@o�@o�@o+@oK�@o\)@o\)@o�P@ol�@o+@p  @p  @o�;@o��@o|�@o+@o+@n�y@nȴ@n��@n�@n�@nȴ@n��@nV@m�T@m/@l�@l�/@m`B@m`B@mp�@m�-@m�-@m�@m/@m�@l�@l��@l�D@l�D@l��@m�@m�@m�@m�@l��@l�j@l�D@l�D@l��@l�j@l��@l�j@l�@l�@l�D@l�@k��@l1@l1@k��@k��@k��@k��@kt�@k@j^5@jn�@j=q@i�#@i�^@i��@i��@i�^@i�^@i��@iX@i��@i��@i��@iX@i&�@i&�@i%@h��@i7L@hĜ@hQ�@hQ�@h �@hb@h  @hb@hb@g�w@g|�@g|�@g|�@g\)@g;d@g;d@g;d@g+@f�y@f�R@f�R@f��@f�+@f�+@f�+@fV@f{@f@e�G�O�@d�@d9X@aG�@`1'@^$�@_
=@`Q�@b�H@d�@f$�@f�R@f�+@dZ@`�u@\1@O�@J�H@J-@:�H@2�!@+t�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBu�Bu�Bu�Bt�Bq�Bo�Bn�Bn�Bp�Bt�Bn�B�B�{B��B��B�'B�LB�RB�XB�qB��BĜBĜBĜBǮBȴBɺB��B��B��B��B��B��B��B��BɺBǮBƨBŢBŢBŢBƨBǮBǮBȴBȴBǮBǮBɺB��BɺBȴBƨBŢBĜBÖBÖBÖBBBBB��B��B��B��B��B��B��B�}B�qB�jB�jB�dB�dB�^B�^B�qB��B��B��B��B��B��B�}B�jB�dB�dB�dB�dB�^B�XB�RB�RB�LB�LB�LB�RB�wB�}B��B��B�}B�}B�wB�qB�qB�wB�wB�wB�wB�wB�wB�qB�qB�jB�jB�dB�dB�dB�dB�dB�dB�^B�^B�XB�XB�XB�RB�LB�LB�XB�XB�XB�RB�RB�RB�RB�RB�RB�RB�RB�LB�LB�LB�LB�RB�XB�^B�dB�jB�jB�dB�dB�^B�XB�RB�XB�XB�XB�RB�LB�LB�LB�FB�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�?B�9B�9B�9B�9B�9B�3B�-B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�!B�-B�3B�3B�3B�3B�3B�9B�9B�3B�3B�3B�3B�3B�-B�!B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�{B�{B�{B�{B�uB�uB�uB�uB�uB�uB�uB�uB�oB�oB�oB�oB�oB�oB�oB�oB�hB�hB�hB�hB�\B�bB�DB�JB�1B�PB�oB��B��B�B�9B�^B�^B�XB�LB��B��B�-B��B��B��111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  BufBuaBueBt[BqKBo@Bn6Bn7BpCBt[Bn5B��B�B�-B��B��B��B��B��B�B�"B�4B�4B�6B�FB�LB�SB�`B�^B�eB�dB�^B�ZB�XB�YB�TB�DB�BB�9B�;B�:B�CB�FB�DB�KB�MB�DB�CB�RB�YB�RB�KB�?B�9B�5B�,B�-B�/B�%B�&B�&B�#B�B�B�B�B�"B�B� B�B�
B�B�B��B��B��B��B�	B�B� B� B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�
B�
B�B�B�B�B�B�B�
B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�zB�oB�oB�rB�{B�xB�zB�rB�{B�sB�zB��B�B�B��B��B��B��B�B�B�{B�xB�{B�xB�{B�yB�yB�zB�tB�sB�rB�sB�qB�rB�|B�B�B��B��B��B��B��B�wB�zB�sB�pB�iB�aB�aB�aB�aB�cB�XB�]B�YB�ZB�UB�ZB�ZB�UB�ZB�[B�YB�ZB�ZB�[B�_B�_B�aB�^B�^B�^B�UB�[B�^B�^B�[B�^B�hB�nB�qB�qB�nB�pB�tB�yB�yB�yB�yB�tB�rB�uB�vB�oB�mB�oB�oB�fB�gB�gB�gB�iB�iB�gB�iB�gB�cB�cB�cB�cB�cB�cB�cB�aB�[B�\B�[B�[B�\B�]B�`B�[B�\B�]B�WB�UB�WB�WB�WB�WB�MB�IB�LB�CB�CB�CB�CB�EB�DB�>B�>B�>B�>B�<B�<B�6B�6B�9B�?B�?B�>B�>B�>B�DB�DB�BB�KB�JB�KB�TB�TB�VB�PB�OB�PB�KB�JB�JB�JB�JB�JB�KB�KB�CB�?B�7B�0B�7B�=B�?B�?B�EB�EB�?B�?B�?B�:B�7B�5B�7B�AB�EB�EB�CB�CB�?B�?B�=B�=B�?B�?B�DB�EB�DB�EB�=B�@B�:B�7B�:B�9B�:B�:B�:B�4B�4B�+B�+B�,B�%B�#B�B�%B�%B�%B�%B�%B�%B�%B�%B�%B�!B�"B�"B�"B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�	B�	B�B� B�B�G�O�B��B��B��B��B��B��B�	B�GB�~B��B��B��B��B��B��B�|B��B��B�VB�EB�D111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt( sw_cndr(PSAL,TEMP,PRES), TEMP, PRES_ADJUSTED )                                                                                                                                                                                         dP =-0.76 dbar.                                                                                                                                                                                                                                                 none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressures adjusted by using pressure offset at the sea surface. The quoted error is manufacturer specified accuracy in dbar.                                                                                                                                    The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   No significant salinity drift detected. The quoted error is max[0.01, 1xOWC uncertainty] in PSS-78.                                                                                                                                                             201904051344042019040513440420190405134404  AO  ARCAADJP                                                                    20181006004031    IP                G�O�G�O�G�O�                AO  ARGQQCPL                                                                    20181006004031  QCP$                G�O�G�O�G�O�5F03E           AO  ARGQQCPL                                                                    20181006004031  QCF$                G�O�G�O�G�O�8000            UW  ARSQUWQC    WOD & nearby Argo as visual check                               20190405134404  IP                  G�O�G�O�G�O�                