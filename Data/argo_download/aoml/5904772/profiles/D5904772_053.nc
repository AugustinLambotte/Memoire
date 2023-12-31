CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS  �   N_CALIB       	N_HISTORY             
   title         Argo float vertical profile    institution       AOML   source        
Argo float     history       2018-10-06T00:40:25Z creation      
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
_FillValue                    �PArgo profile    3.1 1.2 19500101000000  20181006004025  20190405134358  5904772 US ARGO PROJECT                                                 STEPHEN RISER                                                   PRES            TEMP            PSAL               5A   AO  6627                            2C  D   APEX                            7392                            062915                          846 @�/G�g�P1   @�/H`��@N�\)�A�M���1   GPS     Primary sampling: mixed [deeper than nominal 985dbar: discrete; nominal 985dbar to surface: 2dbar-bin averaged]                                                                                                                                                    5A   B   B   @�  @�  A   A   A@  A`  A���A���A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
�fD  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� DhfDh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� DtٚDy�fD�fD�9�D���D�� D�fD�S3D�|�D�� D� D�@ D���D��fD�� D�FfDړ3D��3D�	�D�9�D�y�D��f111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @�  @�  A  A,  AL  Al  A���A���A�  A�  A�  A�  A�  A�  B  B  B  B  B#  B+  B3  B;  BC  BK  BS  B[  Bc  Bk  Bs  B{  B�� B�� B�� B�� B�� B�� B�� B�� B�� B�� B�� B�� B�� B�� B�� B�� B�� Bŀ Bɀ B̀ Bр BՀ Bـ B݀ B� B� B� B� B� B�� B�� B�� C � C� C� C� C� C
� C� C� C� C� C� C� C� C� C� C� C � C"� C$� C&� C(� C*� C,� C.� C0� C2� C4� C6� C8� C:� C<� C>� C@� CB� CD� CF� CH� CJ� CL� CN� CP� CR� CT� CV� CX� CZ� C\� C^� C`� Cb� Cd� Cf� Ch� Cj� Cl� Cn� Cp� Cr� Ct� Cv� Cx� Cz� C|� C~� C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` C�` D 0 D � D0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D	0 D	� D
0 D
�fD0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D0 D� D 0 D � D!0 D!� D"0 D"� D#0 D#� D$0 D$� D%0 D%� D&0 D&� D'0 D'� D(0 D(� D)0 D)� D*0 D*� D+0 D+� D,0 D,� D-0 D-� D.0 D.� D/0 D/� D00 D0� D10 D1� D20 D2� D30 D3� D40 D4� D50 D5� D60 D6� D70 D7� D80 D8� D90 D9� D:0 D:� D;0 D;� D<0 D<� D=0 D=� D>0 D>� D?0 D?� D@0 D@� DA0 DA� DB0 DB� DC0 DC� DD0 DD� DE0 DE� DF0 DF� DG0 DG� DH0 DH� DI0 DI� DJ0 DJ� DK0 DK� DL0 DL� DM0 DM� DN0 DN� DO0 DO� DP0 DP� DQ0 DQ� DR0 DR� DS0 DS� DT0 DT� DU0 DU� DV0 DV� DW0 DW� DX0 DX� DY0 DY� DZ0 DZ� D[0 D[� D\0 D\� D]0 D]� D^0 D^� D_0 D_� D`0 D`� Da0 Da� Db0 Db� Dc0 Dc� Dd0 Dd� De0 De� Df0 Df� Dg0 Dg� Dh6fDh� Di0 Di� Dj0 Dj� Dk0 Dk� Dl0 Dl� Dm0 Dm� Dn0 Dn� Do0 Do� Dp0 Dp� Dq0 Dq� Dr0 Dr� Ds0 Ds� Dt0 Dt� Du	�Dy�fD�fD�Q�D���D�� D�.fD�k3D���D�� D�( D�X D���D��fD� D�^fDګ3D��3D�!�D�Q�D�D��f111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�r�@�r�@�z�@�r�@�r�@�bN@�bN@�Q�@�j@�Q�@�A�@�A�@�A�@�A�@�I�@�I�@�I�@�I�@�I�@�Q�@�Q�@�Q�@�A�@�9X@���@��@ύP@�E�@���@�V@��^@���@�1@�=q@�|�@��@�ff@��h@���@�Q�@��
@��@�$�@�@�V@�Q�@���@�C�@���@�M�@�J@�@��@�G�@��D@��w@��@�l�@�C�@�
=@��@��\@�ff@�-@��-@�7L@���@�Q�@� �@��P@�"�@���@�ff@��@�@��#@��-@���@��h@�p�@�O�@�G�@�x�@���@�%@� �@���@�/@�`B@�p�@�hs@�O�@�?}@�%@���@���@��@��@���@�Ĝ@��j@�j@��m@�K�@�
=@�E�@���@��7@�%@�9X@�  @���@��F@��@�\)@�;d@�@��y@�ȴ@���@�~�@�^5@�=q@�{@�J@���@�p�@�V@��@��/@��/@���@���@���@�Ĝ@��j@��9@��9@��@���@���@���@��@��j@��@��9@��9@��@��@��u@��u@��D@��D@��D@�r�@�r�@�  @|�@|�@�@l�@K�@~��@~�y@~�@~ȴ@~�+@~v�@~ff@~5?@}�T@}��@}O�@|��@|��@|j@|(�@|�@|1@{��@{�m@{�
@{�
@{t�@z�@z�\@z^5@zM�@z=q@z-@zJ@y��@y��@y�7@yX@y7L@y�@x��@x��@x�u@x�u@x�@x�@xr�@xA�@xA�@x �@xb@xb@w�;@w��@w�w@w�P@wl�@wK�@w+@w�@v��@v�@v�R@v�+@vff@vV@v$�@v$�@u�T@u@u��@u�h@u�@t��@t9X@s�m@sƨ@s��@sdZ@s33@s"�@r��@q�@q�@q��@q��@q�7@q�7@q�7@qX@p�u@p  @o��@o+@nV@nE�@n5?@m�@m��@m`B@m?}@m/@mV@l��@mV@mO�@l�/@l��@l�@l��@l��@l�j@l�@l�@l1@k�m@k�m@k�
@k��@k�@k33@k"�@k@j�@j��@j�!@j�!@j��@j�\@j=q@i��@i��@ix�@iX@iG�@iG�@i7L@i7L@i7L@i&�@i&�@i�@i%@h��@h��@h��@h�u@h�@hr�@hbN@hQ�@h1'@h1'@h �@hb@h  @h  @g�@g�@g�@g�@hb@h  @hb@hA�@h  @g��@g�w@g�@g��@g�@g�P@g��@g�;@hb@hA�@hA�@h �@g�@g��@g��@gK�@g�@f�y@f��@f�+@f�+@f�+@fv�@fv�@fV@fV@f$�@fE�@fff@fE�@f{@f@e�@e�@e�@e�@e�@e�T@e��@e��@e��@e��@e��@e��@e��@e��@e��@e��@e��@e��@e�T@f5?@fv�@f��@f��@fff@fE�@e�T@e�T@e�T@e�h@e/@e�@eV@eV@eV@eV@eV@eV@eV@eV@d��@d��@d�@dZ@c�m@c�
@c�m@c�m@dI�@d�j@d�j@d�@d��@dz�@dj@dz�@dj@dj@dZ@dZ@dI�@d(�@c�
@c��@c�@ct�@cdZ@cdZ@cS�@c33@co@c@b�H@b�H@b�@b�@b�@b�@c@b�@b��@b��@c33@co@b�H@b��@b~�@b~�@b~�@b~�@bn�@bn�@b^5@b^5@b=q@a�#@ahs@a�@`bN@_�@_\)@_�@^�+@^��@^ȴ@^�y@^��@^�+@^�+@^ff@^v�@^ff@^v�@^v�@^E�@^@]�-@]�-@]�h@]�@]�@]p�@]O�@]O�@]?}@]/@]/@]�@]V@]V@\��@\��@\��@\��@\��@\�D@[��@Z��@Yhs@XbN@W��@V�y@Vv�@WK�@ZM�@`�@ep�@i&�@hb@h�@e`B@dj@a�@`  @]V@Z-111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  @�r�@�r�@�z�@�r�@�r�@�bN@�bN@�Q�@�j@�Q�@�A�@�A�@�A�@�A�@�I�@�I�@�I�@�I�@�I�@�Q�@�Q�@�Q�@�A�@�9X@���@��@ύP@�E�@���@�V@��^@���@�1@�=q@�|�@��@�ff@��h@���@�Q�@��
@��@�$�@�@�V@�Q�@���@�C�@���@�M�@�J@�@��@�G�@��D@��w@��@�l�@�C�@�
=@��@��\@�ff@�-@��-@�7L@���@�Q�@� �@��P@�"�@���@�ff@��@�@��#@��-@���@��h@�p�@�O�@�G�@�x�@���@�%@� �@���@�/@�`B@�p�@�hs@�O�@�?}@�%@���@���@��@��@���@�Ĝ@��j@�j@��m@�K�@�
=@�E�@���@��7@�%@�9X@�  @���@��F@��@�\)@�;d@�@��y@�ȴ@���@�~�@�^5@�=q@�{@�J@���@�p�@�V@��@��/@��/@���@���@���@�Ĝ@��j@��9@��9@��@���@���@���@��@��j@��@��9@��9@��@��@��u@��u@��D@��D@��D@�r�@�r�@�  @|�@|�@�@l�@K�@~��@~�y@~�@~ȴ@~�+@~v�@~ff@~5?@}�T@}��@}O�@|��@|��@|j@|(�@|�@|1@{��@{�m@{�
@{�
@{t�@z�@z�\@z^5@zM�@z=q@z-@zJ@y��@y��@y�7@yX@y7L@y�@x��@x��@x�u@x�u@x�@x�@xr�@xA�@xA�@x �@xb@xb@w�;@w��@w�w@w�P@wl�@wK�@w+@w�@v��@v�@v�R@v�+@vff@vV@v$�@v$�@u�T@u@u��@u�h@u�@t��@t9X@s�m@sƨ@s��@sdZ@s33@s"�@r��@q�@q�@q��@q��@q�7@q�7@q�7@qX@p�u@p  @o��@o+@nV@nE�@n5?@m�@m��@m`B@m?}@m/@mV@l��@mV@mO�@l�/@l��@l�@l��@l��@l�j@l�@l�@l1@k�m@k�m@k�
@k��@k�@k33@k"�@k@j�@j��@j�!@j�!@j��@j�\@j=q@i��@i��@ix�@iX@iG�@iG�@i7L@i7L@i7L@i&�@i&�@i�@i%@h��@h��@h��@h�u@h�@hr�@hbN@hQ�@h1'@h1'@h �@hb@h  @h  @g�@g�@g�@g�@hb@h  @hb@hA�@h  @g��@g�w@g�@g��@g�@g�P@g��@g�;@hb@hA�@hA�@h �@g�@g��@g��@gK�@g�@f�y@f��@f�+@f�+@f�+@fv�@fv�@fV@fV@f$�@fE�@fff@fE�@f{@f@e�@e�@e�@e�@e�@e�T@e��@e��@e��@e��@e��@e��@e��@e��@e��@e��@e��@e��@e�T@f5?@fv�@f��@f��@fff@fE�@e�T@e�T@e�T@e�h@e/@e�@eV@eV@eV@eV@eV@eV@eV@eV@d��@d��@d�@dZ@c�m@c�
@c�m@c�m@dI�@d�j@d�j@d�@d��@dz�@dj@dz�@dj@dj@dZ@dZ@dI�@d(�@c�
@c��@c�@ct�@cdZ@cdZ@cS�@c33@co@c@b�H@b�H@b�@b�@b�@b�@c@b�@b��@b��@c33@co@b�H@b��@b~�@b~�@b~�@b~�@bn�@bn�@b^5@b^5@b=q@a�#@ahs@a�@`bN@_�@_\)@_�@^�+@^��@^ȴ@^�y@^��@^�+@^�+@^ff@^v�@^ff@^v�@^v�@^E�@^@]�-@]�-@]�h@]�@]�@]p�@]O�@]O�@]?}@]/@]/@]�@]V@]V@\��@\��@\��@\��G�O�@\�D@[��@Z��@Yhs@XbN@W��@V�y@Vv�@WK�@ZM�@`�@ep�@i&�@hb@h�@e`B@dj@a�@`  @]V@Z-111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB
�;B
�;B
�;B
�;B
�5B
�5B
�5B
�5B
�5B
�5B
�5B
�5B
�5B
�5B
�5B
�5B
�5B
�5B
�5B
�5B
�5B
�5B
�/B
�#B
��B
�B�B�B8RBB�BN�BW
BYBZB^5BhsBjBl�Bo�Bp�Bq�Br�Bu�Bv�Bx�B{�B{�B~�B�B�B�B�B�B�B�B�%B�%B�+B�+B�+B�1B�1B�1B�7B�=B�DB�JB�PB�PB�VB�\B�\B�\B�bB�bB�oB�oB�hB�hB�hB�oB�uB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�uB�uB�uB�uB�uB�oB�oB�oB�oB�oB�oB�uB�uB�uB�oB�uB�uB�uB�oB�oB�oB�hB�hB�hB�hB�hB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�\B�\B�\B�\B�\B�\B�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�PB�PB�PB�PB�PB�PB�PB�VB�VB�VB�VB�\B�VB�VB�VB�PB�PB�PB�PB�VB�VB�\B�\B�\B�\B�\B�\B�VB�PB�PB�JB�JB�JB�JB�JB�JB�JB�DB�DB�DB�DB�JB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�JB�JB�PB�PB�JB�JB�JB�DB�DB�DB�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�=B�7B�7B�1B�+B�+B�+B�+B�7B�=B�=B�=B�7B�7B�7B�7B�7B�7B�7B�7B�1B�1B�1B�+B�+B�+B�+B�+B�+B�+B�+B�+B�%B�%B�%B�%B�%B�%B�%B�%B�%B�%B�+B�+B�+B�%B�%B�%B�%B�%B�%B�%B�B�B�B�B�B�B�B� B~�B}�B|�B}�B}�B}�B}�B|�B|�B|�B|�B|�B|�B|�B|�B{�Bz�Bz�Bz�Bz�Bz�Bz�By�By�By�By�By�By�By�By�By�By�By�By�By�Bx�Bv�Bt�Bs�Bp�Bp�Bn�Bm�Bp�By�B�\B��B�B�!B�RB�?B�XB�RB�XB�RB�R111111111111111111111111411111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  B
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
��G�O�B
�-B.B!B7�BB%BNmBV�BX�BY�B]�BhBjBl$Bo6Bp:BqABrEBu[Bv`BxlB{~B{zB~�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B��B�B�B�B�#B�,B�#B�<B�NB�TB�YB�YB�aB�gB�hB�eB�mB�nB�sB�rB�nB�mB�fB�bB�XB�NB�BB�<B�=B�6B�/B�/B�1B�/B�0B�-B�0B�0B�0B�0B�0B�*B�*B�/B�1B�*B�(B�)B�0B�/B�2B�2B�2B�0B�2B�2B�0B�2B�2B�2B�2B�2B�6B�>B�BB�HB�VB�[B�`B�aB�cB�bB�bB�aB�bB�cB�bB�[B�XB�XB�[B�[B�UB�VB�VB�TB�TB�WB�WB�TB�UB�UB�UB�UB�UB�UB�UB�UB�TB�UB�VB�UB�UB�VB�VB�UB�VB�QB�PB�PB�MB�NB�NB�PB�PB�IB�IB�KB�IB�IB�IB�IB�IB�RB�OB�OB�OB�OB�NB�UB�VB�UB�UB�UB�UB�UB�VB�WB�SB�UB�WB�UB�VB�UB�UB�UB�MB�NB�NB�MB�MB�GB�CB�EB�CB�=B�=B�=B�6B�7B�1B�1B�3B�1B�1B�)B�+B�)B�$B�B�B�B�B�B�
B�
B�
B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B� B��B� B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B~�B}�B|�B}�B}�B}�B}�B|�B|�B|�B|�B|�B|�B|�B|�B{�BzzBzxBzyBz{BzyBz{ByuByuByuByuByuByuByuByuByuByuByuByuG�O�BxpBvdBtVBsRBp@Bp@Bn4Bm0Bp@ByyB��B�]B��B��B��B��B��B��B��B��B��111111111111111111111111411111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt( sw_cndr(PSAL,TEMP,PRES), TEMP, PRES_ADJUSTED )                                                                                                                                                                                         dP =-0.75 dbar.                                                                                                                                                                                                                                                 none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressures adjusted by using pressure offset at the sea surface. The quoted error is manufacturer specified accuracy in dbar.                                                                                                                                    The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   No significant salinity drift detected. The quoted error is max[0.01, 1xOWC uncertainty] in PSS-78.                                                                                                                                                             201904051343582019040513435820190405134358  AO  ARCAADJP                                                                    20181006004025    IP                G�O�G�O�G�O�                AO  ARGQQCPL                                                                    20181006004025  QCP$                G�O�G�O�G�O�5F03E           AO  ARGQQCPL                                                                    20181006004025  QCF$                G�O�G�O�G�O�8000            UW  ARSQUWQC    WOD & nearby Argo as visual check                               20190405134358  IP                  G�O�G�O�G�O�                