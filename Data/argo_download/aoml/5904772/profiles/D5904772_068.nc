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
_FillValue                    �PArgo profile    3.1 1.2 19500101000000  20181006004029  20190405134402  5904772 US ARGO PROJECT                                                 STEPHEN RISER                                                   PRES            TEMP            PSAL               DA   AO  6627                            2C  D   APEX                            7392                            062915                          846 @�T��A�1   @�T���� @O���R�@z~��"�1   GPS     Primary sampling: mixed [deeper than nominal 985dbar: discrete; nominal 985dbar to surface: 2dbar-bin averaged]                                                                                                                                                    DA   B   B   @�  @�  A   A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��C�  C�  C�  C�  C��C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'fD'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Dt� Dy��D���D�0 D�� D�� D���D�33D���D�ٚD�  D�,�D���D�ٚD�3D�FfDډ�D��3D�  D�,�D�y�D�|�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @�p�@�p�A
�RA*�RAJ�RAj�RA�\)A�\)A�\)A�\)A�\)A�\)A�\)A�\)B�B
�B�B�B"�B*�B2�B:�BB�BJ�BR�BZ�Bb�Bj�Br�Bz�B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
B�W
C ��C��C��C��C��C
��C��C��C��C��C��C��C��C��C��C��C ��C"��C$��C&��C(��C*��C,��C.��C0��C2��C4��C6��C8��C:��C<��C>��C@��CB��CD��CF��CH��CJ��CL��CN��CP��CR��CT��CV��CX��CZ��C\��C^��C`��Cb��Cd��Cf��Ch��Cj��Cl��Cn��Cp��Cr��Ct��Cv��Cx��Cz��C|��C~��C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�U�C�b�C�U�C�U�C�U�C�U�C�b�C�U�C�U�C�U�C�U�D *�D ��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D	*�D	��D
*�D
��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D*�D��D *�D ��D!*�D!��D"*�D"��D#*�D#��D$*�D$��D%*�D%��D&*�D&��D'1GD'��D(*�D(��D)*�D)��D**�D*��D+*�D+��D,*�D,��D-*�D-��D.*�D.��D/*�D/��D0*�D0��D1*�D1��D2*�D2��D3*�D3��D4*�D4��D5*�D5��D6*�D6��D7*�D7��D8*�D8��D9*�D9��D:*�D:��D;*�D;��D<*�D<��D=*�D=��D>*�D>��D?*�D?��D@*�D@��DA*�DA��DB*�DB��DC*�DC��DD*�DD��DE*�DE��DF*�DF��DG*�DG��DH*�DH��DI*�DI��DJ*�DJ��DK*�DK��DL*�DL��DM*�DM��DN*�DN��DO*�DO��DP*�DP��DQ*�DQ��DR*�DR��DS*�DS��DT*�DT��DU*�DU��DV*�DV��DW*�DW��DX*�DX��DY*�DY��DZ*�DZ��D[*�D[��D\*�D\��D]*�D]��D^*�D^��D_*�D_��D`*�D`��Da*�Da��Db*�Db��Dc*�Dc��Dd*�Dd��De*�De��Df*�Df��Dg*�Dg��Dh*�Dh��Di*�Di��Dj*�Dj��Dk*�Dk��Dl*�Dl��Dm*�Dm��Dn*�Dn��Do*�Do��Dp*�Dp��Dq*�Dq��Dr*�Dr��Ds*�Ds��Dt*�Dt��Du
�Dy��D�D�EqD��qD��qD�D�H�D��D��D�qD�B>D��>D��D��D�[�DڟD��D�qD�B>D�D��>111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�ƨ@�ƨ@��F@��w@�ƨ@��w@��w@��F@�ƨ@�ƨ@�ƨ@�ƨ@�ƨ@��F@���@�l�@�l�@�;d@�v�@�^5@��+@���@�X@�&�@���@��@���@��@��D@�r�@�r�@�j@�A�@�1'@� �@��;@�S�@��T@��7@�/@�Ĝ@�j@��@�;d@��^@�bN@�5?@�G�@�Q�@��@�  @��w@�dZ@�=q@�7L@�9X@��@�ȴ@��-@��@��D@��@�Ĝ@��D@�Q�@��@��;@���@���@�5?@�z�@��9@��\@���@���@�J@�(�@�@��@���@��@��@�l�@�ȴ@��R@���@�ȴ@���@��@���@�V@��j@�j@� �@��@���@��y@���@��@�Ĝ@�I�@��F@�\)@���@���@���@���@��/@���@�j@�;@~��@~��@~V@}�T@}�h@}/@|��@{�m@{dZ@z��@z=q@yx�@w�;@w�@v{@u��@u�-@u�h@u�@u�@up�@up�@uO�@u?}@u?}@u�@tj@s��@so@s@so@sC�@tj@t�/@t�j@t�@t�D@t�D@t9X@s��@s@r�@q��@q�#@q��@q��@q��@q�7@q�7@qx�@qx�@qx�@qhs@qhs@qX@qG�@q7L@q%@p��@p��@p��@pA�@p �@o�;@o�@o�@o�P@o\)@o;d@o;d@o�@o
=@n��@n�y@n��@nV@n5?@n@m@m�@mV@l��@l��@l�@l�D@l�D@l�D@lZ@lI�@l(�@l�@l�@l1@k��@kƨ@k��@k��@k�F@k�m@k�m@k�m@k��@k��@l(�@l9X@l9X@kƨ@k�F@k�m@l�D@l�/@l�@m/@m�h@m��@m�-@m�-@nV@n��@o��@o�@pb@pĜ@p�u@q��@rJ@r�@rJ@q��@q��@q��@q��@q�7@qhs@p��@o�@o��@o��@o��@p�u@p��@p�9@p��@p��@p�`@q%@q%@q%@q%@p��@p�@p�@p�u@p��@p�@pQ�@p1'@o�@o�w@o|�@ol�@o;d@o�@n�+@m�T@m��@m�h@m�h@m�h@m�h@m`B@m/@m�@l�@lz�@lZ@lZ@l�@k��@k�m@kƨ@k��@k�@kdZ@kC�@k33@k33@k33@k"�@ko@k@j�H@j��@j�!@j�!@j~�@jM�@j�@i��@i��@i��@i�7@i��@i��@i��@ix�@iX@iG�@i7L@i7L@i7L@i&�@i�@i&�@i&�@i&�@i�@i%@h��@h��@h�9@h�@hbN@hA�@h �@h  @h  @g�@g��@g�@g�P@g|�@gl�@gl�@gl�@g\)@g\)@gK�@g\)@g\)@g\)@g\)@g\)@g\)@g\)@g\)@gK�@gK�@g;d@g;d@g;d@f��@fȴ@f�R@f��@fv�@fff@fff@fff@fff@fV@fV@f@e�@e��@e�-@e�h@ep�@e`B@eO�@e/@eV@d�@d�@d�/@d�/@d��@d��@d��@d�j@d�@d��@d�D@dz�@dz�@dj@dj@dj@dI�@d9X@d9X@d9X@d(�@d9X@d(�@d1@d1@d1@d1@d�@d�@d�@d(�@d(�@d�@d�@c��@c�m@c��@c��@cS�@c33@co@b�@b�H@b��@b��@b��@b��@b�!@bn�@bn�@bn�@b^5@b^5@b^5@b=q@b�@a�@a�^@ax�@aX@a7L@a�@a%@`��@`��@`�@`r�@`bN@`Q�@`1'@` �@`b@`b@_��@_�@_�w@_�@_��@_��@_�P@_|�@_l�@_\)@_K�@_;d@_;d@_�@_
=@_
=@^�y@^��@^��@^�y@^��@^��@^��@^��@^��@^��@^��@^��@^�y@^�@^�@^�@^�@^ȴ@]�@]V@]�@^v�@[t�@[��@a�7@f�+@i�7@h  @g
=@h1'@g\)@fE�@cƨ@bn�@^v�@]p�@Z^5@W\)111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  @�ƨ@�ƨ@��F@��w@�ƨ@��w@��w@��F@�ƨ@�ƨ@�ƨ@�ƨ@�ƨ@��F@���@�l�@�l�@�;d@�v�@�^5@��+@���@�X@�&�@���@��@���@��@��D@�r�@�r�@�j@�A�@�1'@� �@��;@�S�@��T@��7@�/@�Ĝ@�j@��@�;d@��^@�bN@�5?@�G�@�Q�@��@�  @��w@�dZ@�=q@�7L@�9X@��@�ȴ@��-@��@��D@��@�Ĝ@��D@�Q�@��@��;@���@���@�5?@�z�@��9@��\@���@���@�J@�(�@�@��@���@��@��@�l�@�ȴ@��R@���@�ȴ@���@��@���@�V@��j@�j@� �@��@���@��y@���@��@�Ĝ@�I�@��F@�\)@���@���@���@���@��/@���@�j@�;@~��@~��@~V@}�T@}�h@}/@|��@{�m@{dZ@z��@z=q@yx�@w�;@w�@v{@u��@u�-@u�h@u�@u�@up�@up�@uO�@u?}@u?}@u�@tj@s��@so@s@so@sC�@tj@t�/@t�j@t�@t�D@t�D@t9X@s��@s@r�@q��@q�#@q��@q��@q��@q�7@q�7@qx�@qx�@qx�@qhs@qhs@qX@qG�@q7L@q%@p��@p��@p��@pA�@p �@o�;@o�@o�@o�P@o\)@o;d@o;d@o�@o
=@n��@n�y@n��@nV@n5?@n@m@m�@mV@l��@l��@l�@l�D@l�D@l�D@lZ@lI�@l(�@l�@l�@l1@k��@kƨ@k��@k��@k�F@k�m@k�m@k�m@k��@k��@l(�@l9X@l9X@kƨ@k�F@k�m@l�D@l�/@l�@m/@m�h@m��@m�-@m�-@nV@n��@o��@o�@pb@pĜ@p�u@q��@rJ@r�@rJ@q��@q��@q��@q��@q�7@qhs@p��@o�@o��@o��@o��@p�u@p��@p�9@p��@p��@p�`@q%@q%@q%@q%@p��@p�@p�@p�u@p��@p�@pQ�@p1'@o�@o�w@o|�@ol�@o;d@o�@n�+@m�T@m��@m�h@m�h@m�h@m�h@m`B@m/@m�@l�@lz�@lZ@lZ@l�@k��@k�m@kƨ@k��@k�@kdZ@kC�@k33@k33@k33@k"�@ko@k@j�H@j��@j�!@j�!@j~�@jM�@j�@i��@i��@i��@i�7@i��@i��@i��@ix�@iX@iG�@i7L@i7L@i7L@i&�@i�@i&�@i&�@i&�@i�@i%@h��@h��@h�9@h�@hbN@hA�@h �@h  @h  @g�@g��@g�@g�P@g|�@gl�@gl�@gl�@g\)@g\)@gK�@g\)@g\)@g\)@g\)@g\)@g\)@g\)@g\)@gK�@gK�@g;d@g;d@g;d@f��@fȴ@f�R@f��@fv�@fff@fff@fff@fff@fV@fV@f@e�@e��@e�-@e�h@ep�@e`B@eO�@e/@eV@d�@d�@d�/@d�/@d��@d��@d��@d�j@d�@d��@d�D@dz�@dz�@dj@dj@dj@dI�@d9X@d9X@d9X@d(�@d9X@d(�@d1@d1@d1@d1@d�@d�@d�@d(�@d(�@d�@d�@c��@c�m@c��@c��@cS�@c33@co@b�@b�H@b��@b��@b��@b��@b�!@bn�@bn�@bn�@b^5@b^5@b^5@b=q@b�@a�@a�^@ax�@aX@a7L@a�@a%@`��@`��@`�@`r�@`bN@`Q�@`1'@` �@`b@`b@_��@_�@_�w@_�@_��@_��@_�P@_|�@_l�@_\)@_K�@_;d@_;d@_�@_
=@_
=@^�y@^��@^��@^�y@^��@^��@^��@^��@^��@^��@^��@^��@^�y@^�@^�@^�G�O�@^ȴ@]�@]V@]�@^v�@[t�@[��@a�7@f�+@i�7@h  @g
=@h1'@g\)@fE�@cƨ@bn�@^v�@]p�@Z^5@W\)111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�qB�jB�jB�jB�jB�jB�jB�dB�dB�dB�dB�dB�dB�dB�^B�^B�^B�^B�^B�^B�^B�XB�XB�LB�FB�?B�9B�9B�3B�-B�!B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B��B��B��B��B��B��B��B��B�B�'B�'B�B��B��B��B��B��B��B�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�{B�{B�uB�uB�oB�oB�hB�bB�bB�bB�bB�bB�bB�bB�bB�hB�hB�hB�hB�bB�\B�VB�VB�VB�VB�\B�oB�uB�uB�oB�oB�oB�hB�bB�\B�\B�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�PB�PB�PB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�DB�DB�DB�DB�=B�=B�7B�7B�7B�7B�7B�=B�7B�7B�7B�7B�7B�7B�7B�7B�1B�1B�7B�7B�7B�=B�=B�=B�=B�=B�DB�=B�=B�DB�JB�PB�PB�VB�\B�\B�\B�\B�hB�uB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�{B�uB�{B�uB�uB�uB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�hB�hB�hB�hB�hB�hB�hB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�VB�VB�VB�VB�VB�PB�PB�PB�PB�PB�PB�PB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�JB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�JB�JB�JB�JB�JB�JB�DB�DB�DB�DB�=B�=B�=B�=B�=B�7B�7B�7B�7B�7B�7B�7B�7B�7B�7B�1B�1B�1B�+B�+B�+B�%B�%B�%B�%B�%B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B� B� B�B�B~�B� B�bB��B�B��B�B�LB�^B�jB�jB�qB�jB�qB�jB�j111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B� B� B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�~B�B��B��B��B��B��B��B��B��B��B��B�tB�~B��B��B��B��B��B�oB�SB�<B�*B�"B� B�)B�/B�)B�)B�/B�4B�:B�AB�AB�AB�AB�JB�IB�HB�@B�BB�BB�;B�;B�;B�4B�4B�/B�/B�.B�)B�&B�&B�"B�#B�"B�!B�B�B�B�B�B�B�B�B�
B�B�B�B�B�B�B�B�B�	B�	B�	B�	B�B��B��B��B��B��B��B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�
B�B�B�%B�%B�)B�(B�<B�<B�<B�@B�CB�CB�@B�@B�>B�:B�6B�-B�)B�)B�.B�;B�9B�<B�CB�EB�@B�@B�KB�IB�@B�CB�CB�IB�HB�IB�IB�KB�HB�@B�CB�=B�=B�=B�4B�5B�.B�.B�-B�.B�-B�.B�3B�/B�*B�*B�+B�%B�%B�"B�"B�$B�B�B�B�B�B�B�B�B�B�!B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B�B��B�B�B��B��B��B��B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��G�O�B��B�B�B��B��B~�B�B�
B�nB��B��B��B��B�B�B�B�B�B�B�B�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt( sw_cndr(PSAL,TEMP,PRES), TEMP, PRES_ADJUSTED )                                                                                                                                                                                         dP =-0.67 dbar.                                                                                                                                                                                                                                                 none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressures adjusted by using pressure offset at the sea surface. The quoted error is manufacturer specified accuracy in dbar.                                                                                                                                    The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   No significant salinity drift detected. The quoted error is max[0.01, 1xOWC uncertainty] in PSS-78.                                                                                                                                                             201904051344022019040513440220190405134402  AO  ARCAADJP                                                                    20181006004029    IP                G�O�G�O�G�O�                AO  ARGQQCPL                                                                    20181006004029  QCP$                G�O�G�O�G�O�5F03E           AO  ARGQQCPL                                                                    20181006004029  QCF$                G�O�G�O�G�O�8000            UW  ARSQUWQC    WOD & nearby Argo as visual check                               20190405134402  IP                  G�O�G�O�G�O�                