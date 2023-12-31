CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   N_CALIB       	N_HISTORY             
   title         Argo float vertical profile    institution       AOML   source        
Argo float     history       2018-10-06T00:40:28Z creation      
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
_FillValue                    �tArgo profile    3.1 1.2 19500101000000  20181006004028  20190405134401  5904772 US ARGO PROJECT                                                 STEPHEN RISER                                                   PRES            TEMP            PSAL               AA   AO  6627                            2C  D   APEX                            7392                            062915                          846 @�M;Z���1   @�M<��@N�ȴ9X�@���v�1   GPS     Primary sampling: mixed [deeper than nominal 985dbar: discrete; nominal 985dbar to surface: 2dbar-bin averaged]                                                                                                                                                    AA   B   B   @&ff@�  @�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bg��Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C�C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� DgfDg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Dt�3Dy�3D���D�@ D�|�D���D�	�D�@ D�ffD�ɚD� D�9�D�i�D��fD��D�9�D�|�D�� D�3D�6fD�l�D��31111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @XQ�@���@���Az�A,z�ALz�Alz�A�=qA�=qA�=qA�=qA�=qA�=qA�=qA�=qB�B�B�B�B#�B+�B3�B;�BC�BK�BS�B[�Bc�Bj�RBs�B{�B��\B��\B��\B��\B��\B��\B��\B��\B��\B��\B��\B��\B�\)B��\B��\B��\B��\Bŏ\Bɏ\B͏\Bя\B�\)Bُ\Bݏ\B�\B�\B�\B�\B�\B��\B��\B��\C ǮCǮCǮCǮC�HC
ǮCǮCǮCǮCǮCǮCǮCǮCǮCǮCǮC ǮC"ǮC$ǮC&ǮC(ǮC*ǮC,ǮC.ǮC0ǮC2ǮC4ǮC6ǮC8ǮC:ǮC<ǮC>ǮC@ǮCBǮCDǮCFǮCHǮCJǮCLǮCNǮCPǮCRǮCTǮCVǮCXǮCZǮC\ǮC^ǮC`ǮCbǮCdǮCfǮChǮCjǮClǮCnǮCpǮCrǮCtǮCvǮCxǮCzǮC|ǮC~ǮC�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�p�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�C�c�D 1�D ��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D	1�D	��D
1�D
��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D1�D��D 1�D ��D!1�D!��D"1�D"��D#1�D#��D$1�D$��D%1�D%��D&1�D&��D'1�D'��D(1�D(��D)1�D)��D*1�D*��D+1�D+��D,1�D,��D-1�D-��D.1�D.��D/1�D/��D01�D0��D11�D1��D21�D2��D31�D3��D41�D4��D51�D5��D61�D6��D71�D7��D81�D8��D91�D9��D:1�D:��D;1�D;��D<1�D<��D=1�D=��D>1�D>��D?1�D?��D@1�D@��DA1�DA��DB1�DB��DC1�DC��DD1�DD��DE1�DE��DF1�DF��DG1�DG��DH1�DH��DI1�DI��DJ1�DJ��DK1�DK��DL1�DL��DM1�DM��DN1�DN��DO1�DO��DP1�DP��DQ1�DQ��DR1�DR��DS1�DS��DT1�DT��DU1�DU��DV1�DV��DW1�DW��DX1�DX��DY1�DY��DZ1�DZ��D[1�D[��D\1�D\��D]1�D]��D^1�D^��D_1�D_��D`1�D`��Da1�Da��Db1�Db��Dc1�Dc��Dd1�Dd��De1�De��Df1�Df��Dg8RDg��Dh1�Dh��Di1�Di��Dj1�Dj��Dk1�Dk��Dl1�Dl��Dm1�Dm��Dn1�Dn��Do1�Do��Dp1�Dp��Dq1�Dq��Dr1�Dr��Ds1�Ds��Dt1�Dt��DuDy�D��D�X�D���D���D�"�D�X�D�\D��D�(�D�R�D���D��\D�%�D�R�Dڕ�D���D�,)D�O\D��D��)1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@\)@\)@\)@\)@l�@l�@l�@|�@|�@�P@|�@�P@��@�P@��@�@��@��@��@�P@��@�P@�P@��@�@�@�@�@�@�@�@�@�w@�w@�@�@��@��@�w@�w@�w@��@��@��@�;@�;@�;@�;@�;@�@�@�@�@�@�@�  @�  @�  @�  @�  @�1@�1@�1@�1@�1@�b@�b@�b@�b@��@��@��@��@��@��@��@��@� �@� �@��@�(�@�1'@� �@� �@�(�@�9X@�1'@�1'@�1'@�9X@�9X@�9X@�9X@�9X@�1'@�1'@�(�@�(�@�(�@�(�@�1'@�1'@�1'@�9X@�Q�@�Q�@�Z@�Z@�Z@�Z@�Z@�Z@�bN@�Z@�bN@�bN@�bN@�bN@�bN@�bN@�bN@�bN@�bN@�bN@�bN@�bN@�j@�j@�j@�j@�j@�j@�j@�j@�r�@�r�@�r�@�r�@�z�@�r�@�z�@�z�@��@��@��@��@��@��@��@��@��D@��D@��D@��D@��D@��D@��D@��u@��u@��u@��u@��u@��u@���@���@���@���@���@���@��@��@��@��9@��9@��9@��9@��9@��j@��j@��j@��9@���@���@��@��9@��j@��j@���@��D@�z�@�r�@�j@�r�@�z�@�r�@�bN@�j@�r�@�j@�1'@� �@�;@~�+@}�@}��@}�-@}��@}`B@|�@|z�@|j@|Z@|�@{�m@{�F@{�@{��@{"�@z�!@z^5@z=q@z�@z�@zJ@y��@y��@y�@y�#@y�#@y�#@y��@y��@y�#@|1@|��@}�@}`B@}�T@}�@}�@}�@}�@}�@~{@~$�@~ff@�@l�@\)@�b@�1'@�9X@K�@~��@~ff@~5?@~5?@~E�@~5?@~$�@~@}`B@|��@|�@|��@{�F@{S�@zM�@yhs@xĜ@x�u@w�w@w�P@wl�@v�y@v5?@v@v@u�@u�T@u��@u��@u@u@u@u��@u�@tj@sƨ@s��@s��@sdZ@s@r��@r-@q�^@qX@qX@qX@qG�@q7L@q&�@q%@p�@pbN@p �@pb@pb@p  @o�w@o�P@o�P@o�P@o�P@ol�@ol�@ol�@oK�@o;d@o;d@o+@o�@o
=@n�@nȴ@nȴ@n��@nv�@nV@n5?@n5?@n@m�@m��@m@m�-@m�h@mp�@mO�@m�@l�/@l��@lj@l1@k�m@k�m@k�m@k�m@k�@ko@k@k@k@k@j�@j�H@j�H@j��@j�!@j�!@j��@j��@j��@j�\@j�\@j�\@j~�@j~�@j~�@jn�@jn�@jn�@jn�@jn�@j~�@j~�@j~�@j�\@j�\@j��@j��@j�!@j��@j��@j��@j��@j�!@j��@j��@j�!@jn�@j=q@j�@i�#@i�#@i��@i��@i�^@i�^@i�^@i��@i��@i��@i��@i�^@i�7@i7L@iX@iX@iG�@i&�@i�@i�@i�@h��@h��@h��@h�@hr�@hbN@hA�@hA�@h1'@h �@h �@h �@h �@h �@hb@hb@hb@g�@g�w@g��@g|�@g\)@g\)@g;d@g+@g�@g
=@g
=@f��@f��@f�y@f�R@f�+@f�+@fE�@fE�@fE�@f5?@f5?@f5?@f$�@e�T@e�@ep�@e`B@eO�@eO�@eO�@e/@e�@eV@d�@d�/@d��@d�D@d�D@dz�@d�D@d��@d�j@d��@eV@d�/@d��@d�@d�@d�@d�/@d�/@d��@d�j@d��@d��@d��@d��@d��@dj@b�!@`r�@_�@^��@\�@\1@[t�@[@]?}@a�@f��@g��@i�^@iX@g�;@ep�@b�H@`  @]�h@[331111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 @\)@\)@\)@\)@l�@l�@l�@|�@|�@�P@|�@�P@��@�P@��@�@��@��@��@�P@��@�P@�P@��@�@�@�@�@�@�@�@�@�w@�w@�@�@��@��@�w@�w@�w@��@��@��@�;@�;@�;@�;@�;@�@�@�@�@�@�@�  @�  @�  @�  @�  @�1@�1@�1@�1@�1@�b@�b@�b@�b@��@��@��@��@��@��@��@��@� �@� �@��@�(�@�1'@� �@� �@�(�@�9X@�1'@�1'@�1'@�9X@�9X@�9X@�9X@�9X@�1'@�1'@�(�@�(�@�(�@�(�@�1'@�1'@�1'@�9X@�Q�@�Q�@�Z@�Z@�Z@�Z@�Z@�Z@�bN@�Z@�bN@�bN@�bN@�bN@�bN@�bN@�bN@�bN@�bN@�bN@�bN@�bN@�j@�j@�j@�j@�j@�j@�j@�j@�r�@�r�@�r�@�r�@�z�@�r�@�z�@�z�@��@��@��@��@��@��@��@��@��D@��D@��D@��D@��D@��D@��D@��u@��u@��u@��u@��u@��u@���@���@���@���@���@���@��@��@��@��9@��9@��9@��9@��9@��j@��j@��j@��9@���@���@��@��9@��j@��j@���@��D@�z�@�r�@�j@�r�@�z�@�r�@�bN@�j@�r�@�j@�1'@� �@�;@~�+@}�@}��@}�-@}��@}`B@|�@|z�@|j@|Z@|�@{�m@{�F@{�@{��@{"�@z�!@z^5@z=q@z�@z�@zJ@y��@y��@y�@y�#@y�#@y�#@y��@y��@y�#@|1@|��@}�@}`B@}�T@}�@}�@}�@}�@}�@~{@~$�@~ff@�@l�@\)@�b@�1'@�9X@K�@~��@~ff@~5?@~5?@~E�@~5?@~$�@~@}`B@|��@|�@|��@{�F@{S�@zM�@yhs@xĜ@x�u@w�w@w�P@wl�@v�y@v5?@v@v@u�@u�T@u��@u��@u@u@u@u��@u�@tj@sƨ@s��@s��@sdZ@s@r��@r-@q�^@qX@qX@qX@qG�@q7L@q&�@q%@p�@pbN@p �@pb@pb@p  @o�w@o�P@o�P@o�P@o�P@ol�@ol�@ol�@oK�@o;d@o;d@o+@o�@o
=@n�@nȴ@nȴ@n��@nv�@nV@n5?@n5?@n@m�@m��@m@m�-@m�h@mp�@mO�@m�@l�/@l��@lj@l1@k�m@k�m@k�m@k�m@k�@ko@k@k@k@k@j�@j�H@j�H@j��@j�!@j�!@j��@j��@j��@j�\@j�\@j�\@j~�@j~�@j~�@jn�@jn�@jn�@jn�@jn�@j~�@j~�@j~�@j�\@j�\@j��@j��@j�!@j��@j��@j��@j��@j�!@j��@j��@j�!@jn�@j=q@j�@i�#@i�#@i��@i��@i�^@i�^@i�^@i��@i��@i��@i��@i�^@i�7@i7L@iX@iX@iG�@i&�@i�@i�@i�@h��@h��@h��@h�@hr�@hbN@hA�@hA�@h1'@h �@h �@h �@h �@h �@hb@hb@hb@g�@g�w@g��@g|�@g\)@g\)@g;d@g+@g�@g
=@g
=@f��@f��@f�y@f�R@f�+@f�+@fE�@fE�@fE�@f5?@f5?@f5?@f$�@e�T@e�@ep�@e`B@eO�@eO�@eO�@e/@e�@eV@d�@d�/@d��@d�D@d�D@dz�@d�D@d��@d�j@d��@eV@d�/@d��@d�@d�@d�@d�/@d�/@d��@d�j@d��@d��@d��@d��G�O�@dj@b�!@`r�@_�@^��@\�@\1@[t�@[@]?}@a�@f��@g��@i�^@iX@g�;@ep�@b�H@`  @]�h@[331111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�=B�DB�=B�DB�=B�DB�DB�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�7B�7B�=B�=B�7B�7B�7B�7B�7B�7B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�DB�=B�=B�=B�=B�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�=B�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�JB�JB�JB�PB�VB�\B�\B�\B�VB�VB�bB�hB�oB�oB�oB�oB�oB�{B��B��B��B��B��B��B��B�{B�uB�uB�uB�uB�oB�oB�hB�hB�hB�bB�bB�bB�bB�\B�\B�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�\B�bB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�{B�{B�{B�{B�uB�uB�oB�uB�uB�oB�oB�uB�uB�uB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�uB�uB�uB�uB�uB�uB�{B�{B�{B�{B�{B�{B�{B�uB�uB�uB�oB�uB�uB�oB�oB�oB�oB�oB�oB�oB�uB�uB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�bB�bB�bB�bB�bB�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�VB�VB�VB�PB�PB�PB�PB�PB�PB�PB�JB�JB�JB�JB�JB�DB�DB�DB�DB�DB�DB�DB�DB�=B�=B�=B�=B�=B�DB�DB�DB�DB�DB�DB�JB�JB�JB�DB�DB�DB�DB�DB�DB�DB�DB�DB�=B�1B�B�B�B~�B}�B|�B{�B�B�hB��B�B�FB�LB�dB�dB�dB�jB�dB�j1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�	B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�'B�.B�;B�;B�;B�?B�?B�EB�KB�SB�WB�dB�fB�nB�vB�}B��B�~B�yB�vB�vB�uB�vB�wB�vB�vB�pB�pB�sB�pB�kB�dB�_B�RB�PB�UB�LB�KB�NB�EB�FB�GB�FB�?B�@B�=B�?B�>B�@B�=B�=B�:B�<B�4B�2B�2B�2B�4B�.B�,B�,B�'B�(B�'B�(B�(B�'B�'B�(B�(B�#B�!B�!B�!B�!B�!B�!B� B�!B� B�B�!B� B�!B� B� B�#B�!B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�
B�B�B�B�
B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�
B�
B�
B�B�B�B�B�B�B�B�
B�B�B�B�B�B�B�B�B�B�B�B�B�B�	B�B�B�B�B�B�B�B�B�B�B�B��B�B� B� B� B� B��B� B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��G�O�B��B��B��B��B��B~�B}�B|�B{�B��B�B�wB��B��B��B� B�B�B�B�B�1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111 <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt( sw_cndr(PSAL,TEMP,PRES), TEMP, PRES_ADJUSTED )                                                                                                                                                                                         dP =-0.78 dbar.                                                                                                                                                                                                                                                 none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressures adjusted by using pressure offset at the sea surface. The quoted error is manufacturer specified accuracy in dbar.                                                                                                                                    The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   No significant salinity drift detected. The quoted error is max[0.01, 1xOWC uncertainty] in PSS-78.                                                                                                                                                             201904051344012019040513440120190405134401  AO  ARCAADJP                                                                    20181006004028    IP                G�O�G�O�G�O�                AO  ARGQQCPL                                                                    20181006004028  QCP$                G�O�G�O�G�O�5F03E           AO  ARGQQCPL                                                                    20181006004028  QCF$                G�O�G�O�G�O�8000            UW  ARSQUWQC    WOD & nearby Argo as visual check                               20190405134401  IP                  G�O�G�O�G�O�                