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
_FillValue                    �PArgo profile    3.1 1.2 19500101000000  20181006004028  20190405134401  5904772 US ARGO PROJECT                                                 STEPHEN RISER                                                   PRES            TEMP            PSAL               ?A   AO  6627                            2C  D   APEX                            7392                            062915                          846 @�H=�y�1   @�H>m�T2@N˅�Q��A�dZ�1   GPS     Primary sampling: mixed [deeper than nominal 985dbar: discrete; nominal 985dbar to surface: 2dbar-bin averaged]                                                                                                                                                    ?A   B   B   @�  @�  A   A   AA��A`  A�  A�  A�  A�  A���A�  A�  A�  B   B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B`  Bh  Bp  Bx  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB�CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CZ  C\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D#��D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/  D/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DV� DW  DW� DX  DX� DY  DY� DZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� DtٚDy�3D��D�@ D�|�D�ٚD�  D�C3D�y�D��3D�  D�0 D��fD�� D� D�,�Dڃ3D��3D���D�I�D�p D���111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @�\)@�\)A�A+�AMG�Ak�A��
A��
A��
A��
Aƣ�A��
A��
A��
B�B
�B�B�B"�B*�B2�B:�BB�BJ�BR�BZ�Bb�Bj�Br�Bz�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�C ��C��C��C��C��C
��C��C��C��C��C��C��C��C��C��C��C ��C"��C$��C&��C(��C*��C,��C.��C0��C2��C4��C6��C8��C:��C<��C>��C@��CB�{CD��CF��CH��CJ��CL��CN��CP��CR��CT��CV��CX��CZ��C\��C^��C`��Cb��Cd��Cf��Ch��Cj��Cl��Cn��Cp��Cr��Ct��Cv��Cx��Cz��C|��C~��C�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qD .�D ��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D	.�D	��D
.�D
��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D .�D ��D!.�D!��D".�D"��D#.�D#��D$(RD$��D%.�D%��D&.�D&��D'.�D'��D(.�D(��D).�D)��D*.�D*��D+.�D+��D,.�D,��D-.�D-��D..�D.��D/.�D/��D0.�D0��D1.�D1��D2.�D2��D3.�D3��D4.�D4��D5.�D5��D6.�D6��D7.�D7��D8.�D8��D9.�D9��D:.�D:��D;.�D;��D<.�D<��D=.�D=��D>.�D>��D?.�D?��D@.�D@��DA.�DA��DB.�DB��DC.�DC��DD.�DD��DE.�DE��DF.�DF��DG.�DG��DH.�DH��DI.�DI��DJ.�DJ��DK.�DK��DL.�DL��DM.�DM��DN.�DN��DO.�DO��DP.�DP��DQ.�DQ��DR.�DR��DS.�DS��DT.�DT��DU.�DU��DV.�DV��DW.�DW��DX.�DX��DY.�DY��DZ.�DZ��D[.�D[��D\.�D\��D].�D]��D^.�D^��D_.�D_��D`.�D`��Da.�Da��Db.�Db��Dc.�Dc��Dd.�Dd��De.�De��Df.�Df��Dg.�Dg��Dh.�Dh��Di.�Di��Dj.�Dj��Dk.�Dk��Dl.�Dl��Dm.�Dm��Dn.�Dn��Do.�Do��Dp.�Dp��Dq.�Dq��Dr.�Dr��Ds.�Ds��Dt.�Dt��DuRDy��D�$)D�W\D��)D���D�\D�Z�D���D�ڏD�\D�G\D���D��\D�'\D�D)Dښ�D�ڏD�)D�`�D�\D���111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@���@�@�@�@�@���@���@���@���@���@���@���@��#@��#@��@��T@��@��@��#@��T@��#@���@��#@���@���@���@�@�@�@�@�@�@�@�J@�J@�J@�{@�{@�{@�{@�{@�J@�J@�{@�{@�{@�{@�{@�{@��@��@�$�@�$�@�-@�$�@�$�@�$�@�$�@�-@�-@�-@�5?@�5?@�5?@�5?@�=q@�5?@�5?@�5?@�$�@�$�@�-@�-@�-@�-@�-@�-@�-@�-@�5?@�=q@�5?@�=q@�=q@�=q@�E�@�=q@�M�@�5?@�{@��@��@��@��@��#@��#@���@���@���@��#@��#@��#@��#@��#@��#@��T@��T@��T@��@�J@�v�@��R@��y@�+@��@��@�V@��@�hs@�X@�`B@�hs@�x�@�x�@��@��7@���@�-@���@�ȴ@�ȴ@��+@�;d@�33@�n�@�V@�M�@�5?@�=q@�E�@��@��@�n�@��y@�;d@�+@�
=@��@��@�ȴ@��\@�$�@��^@��@���@�9X@�b@���@���@��@�o@�ȴ@��+@���@���@�n�@�=q@�@��^@�&�@��@�%@��u@�j@�1'@��
@�t�@���@�ȴ@�~�@�n�@��@��h@��h@�/@��`@�r�@�(�@��m@��@���@�t�@�C�@��@��y@���@��R@���@��\@��+@�^5@��@��#@���@��@�x�@�x�@�x�@�hs@�/@�V@��`@���@�Ĝ@��j@��@��@�r�@�9X@��@�b@���@���@��@���@�dZ@�S�@�S�@�K�@��@��R@�~�@�^5@�V@�E�@�5?@�$�@�{@�J@���@���@�@��h@�p�@�/@�%@�Ĝ@��u@�A�@�@��@+@}�h@}�@|�@|z�@|I�@|(�@{�m@{�
@{�
@{��@{�@{33@z�H@z�!@zM�@z�@zJ@y��@y7L@x�u@x �@w�w@wl�@v�y@vff@v5?@u�@u@u��@u`B@u�@t�j@tj@tI�@t(�@s��@sƨ@sdZ@r�@r�!@rn�@rM�@r=q@r-@r�@q��@q�@q�@q�#@q�#@q��@q�^@q��@q�7@qG�@q&�@q�@p�`@pĜ@p��@p�@pr�@pbN@pA�@p1'@p �@pb@p  @p  @p  @o�@o�;@o��@o�w@o�w@o��@o�;@o�;@o�;@o�w@o��@o|�@o\)@o;d@o;d@o;d@o;d@o;d@o+@o�@o
=@n�y@n�R@n�+@nff@nE�@n$�@n5?@nE�@n5?@n{@n{@m�@m��@m@m@m@m@m�-@m@m@m@m��@m��@m��@m��@m�@mp�@mO�@m�@mV@l��@l��@l�@l��@lz�@lZ@lI�@l(�@l�@l1@k�m@k�
@k�F@k�F@kƨ@kƨ@k�F@k��@k��@kt�@kdZ@kdZ@kdZ@kS�@kC�@ko@k@j�@k@j��@j~�@j^5@j=q@j=q@j=q@j=q@j-@j�@jJ@i�#@i��@i�^@i�7@i�7@ix�@ix�@ix�@ihs@ihs@iG�@i7L@i�@i%@i�@h�`@hĜ@h�9@h�9@h�9@h��@h��@h��@h�u@hr�@h1'@h �@h �@h �@g�@g�@g�;@g�;@g��@g��@g��@g��@g�w@g�@g�P@g|�@g\)@gK�@gK�@gK�@gK�@g;d@g+@g+@g+@g+@g+@g�@f��@f�y@fȴ@f�R@f�R@f�R@f�R@f�R@f��@f��@f�+@fv�@fff@fE�@f$�@f{@f{@f@f@f@e�T@e�-@e@e@d�@c��@a�^@`bN@_
=@]p�@[@Z��@[ƨ@[C�@]`B@a��@f5?@g�@i��@i�7@g��@e?}@b�H@`Ĝ@^111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  @���@�@�@�@�@���@���@���@���@���@���@���@��#@��#@��@��T@��@��@��#@��T@��#@���@��#@���@���@���@�@�@�@�@�@�@�@�J@�J@�J@�{@�{@�{@�{@�{@�J@�J@�{@�{@�{@�{@�{@�{@��@��@�$�@�$�@�-@�$�@�$�@�$�@�$�@�-@�-@�-@�5?@�5?@�5?@�5?@�=q@�5?@�5?@�5?@�$�@�$�@�-@�-@�-@�-@�-@�-@�-@�-@�5?@�=q@�5?@�=q@�=q@�=q@�E�@�=q@�M�@�5?@�{@��@��@��@��@��#@��#@���@���@���@��#@��#@��#@��#@��#@��#@��T@��T@��T@��@�J@�v�@��R@��y@�+@��@��@�V@��@�hs@�X@�`B@�hs@�x�@�x�@��@��7@���@�-@���@�ȴ@�ȴ@��+@�;d@�33@�n�@�V@�M�@�5?@�=q@�E�@��@��@�n�@��y@�;d@�+@�
=@��@��@�ȴ@��\@�$�@��^@��@���@�9X@�b@���@���@��@�o@�ȴ@��+@���@���@�n�@�=q@�@��^@�&�@��@�%@��u@�j@�1'@��
@�t�@���@�ȴ@�~�@�n�@��@��h@��h@�/@��`@�r�@�(�@��m@��@���@�t�@�C�@��@��y@���@��R@���@��\@��+@�^5@��@��#@���@��@�x�@�x�@�x�@�hs@�/@�V@��`@���@�Ĝ@��j@��@��@�r�@�9X@��@�b@���@���@��@���@�dZ@�S�@�S�@�K�@��@��R@�~�@�^5@�V@�E�@�5?@�$�@�{@�J@���@���@�@��h@�p�@�/@�%@�Ĝ@��u@�A�@�@��@+@}�h@}�@|�@|z�@|I�@|(�@{�m@{�
@{�
@{��@{�@{33@z�H@z�!@zM�@z�@zJ@y��@y7L@x�u@x �@w�w@wl�@v�y@vff@v5?@u�@u@u��@u`B@u�@t�j@tj@tI�@t(�@s��@sƨ@sdZ@r�@r�!@rn�@rM�@r=q@r-@r�@q��@q�@q�@q�#@q�#@q��@q�^@q��@q�7@qG�@q&�@q�@p�`@pĜ@p��@p�@pr�@pbN@pA�@p1'@p �@pb@p  @p  @p  @o�@o�;@o��@o�w@o�w@o��@o�;@o�;@o�;@o�w@o��@o|�@o\)@o;d@o;d@o;d@o;d@o;d@o+@o�@o
=@n�y@n�R@n�+@nff@nE�@n$�@n5?@nE�@n5?@n{@n{@m�@m��@m@m@m@m@m�-@m@m@m@m��@m��@m��@m��@m�@mp�@mO�@m�@mV@l��@l��@l�@l��@lz�@lZ@lI�@l(�@l�@l1@k�m@k�
@k�F@k�F@kƨ@kƨ@k�F@k��@k��@kt�@kdZ@kdZ@kdZ@kS�@kC�@ko@k@j�@k@j��@j~�@j^5@j=q@j=q@j=q@j=q@j-@j�@jJ@i�#@i��@i�^@i�7@i�7@ix�@ix�@ix�@ihs@ihs@iG�@i7L@i�@i%@i�@h�`@hĜ@h�9@h�9@h�9@h��@h��@h��@h�u@hr�@h1'@h �@h �@h �@g�@g�@g�;@g�;@g��@g��@g��@g��@g�w@g�@g�P@g|�@g\)@gK�@gK�@gK�@gK�@g;d@g+@g+@g+@g+@g+@g�@f��@f�y@fȴ@f�R@f�R@f�R@f�R@f�R@f��@f��@f�+@fv�@fff@fE�@f$�@f{@f{@f@f@f@e�T@e�-@eG�O�@d�@c��@a�^@`bN@_
=@]p�@[@Z��@[ƨ@[C�@]`B@a��@f5?@g�@i��@i�7@g��@e?}@b�H@`Ĝ@^111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oBq�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bp�Bq�Bq�Bq�Bq�Bq�Bq�Br�Br�Bq�Bq�Br�Br�Br�Br�Br�Br�Br�Br�Br�Br�Bs�Bs�Bw�B{�B� B�B�JB�bB�uB��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�!B�-B�?B�XB�qB�wB�wB�wB�wB�wB�qB�qB�qB�jB�jB�dB�dB�dB�^B�^B�XB�RB�RB�RB�RB�RB�LB�LB�FB�?B�?B�9B�9B�3B�3B�-B�-B�'B�!B�!B�!B�!B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��BdZB|�B��BW
B�#BcTB�PB~�B�B~�Bu�BjBaHBdZBjBl�Bn�Bp�Bt�By�Bz�B{�B|�B{�B{�B|�B~�B�B� B~�B}�B|�B� B�%B�1B�=B�DB�DB�DB�DB�DB�DB�=B�=B�1B�+B�%B�%B�B�+B�1B�1B�7B�7B�=B�DB�DB�DB�JB�PB�PB�VB�PB�PB�PB�PB�PB�VB�VB�VB�PB�=B�1B�+B�1B�1B�=B�=B�JB�DB�=B�7B�1B�+B�+B�%B�+B�7B�JB�PB�\B�oB�bB�\B�PB�JB�JB�PB�\B�bB�hB�hB�hB�oB�oB�uB�{B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�{B�{B�{B�{B�{B�{B�uB�uB�uB�uB�uB�uB�oB�oB�oB�oB�oB�oB�oB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�bB�bB�bB�bB�bB�bB�bB�bB�bB�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�VB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�JB�JB�JB�JB�JB�JB�JB�JB�DB�DB�DB�DB�=B�+B�B�B�B|�Bw�Bx�B{�Bz�B�B�bB��B��B�9B�XB�dB�^B�qB�qB�q411111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111144444444111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  G�O�Bp=Bp<Bp<Bp=Bp=Bp;Bp=Bp;Bp?Bp<Bp<Bp<Bp?Bp<Bp?Bp;Bp<Bp=Bp=Bp?Bp=Bp?Bp?Bp=Bp=Bp;Bp;Bp;Bp=Bp<Bp<Bp;Bp?Bp=Bp=Bp<Bp>Bp=Bp?Bp<Bp<Bp<Bp>Bp=Bp>Bp>Bp>Bp=Bp>Bp>Bp;Bp<Bp<Bp;Bp<Bp;Bp<Bp?Bp;Bp>Bp<Bp<Bp<Bp<Bp<Bp>Bp<Bp<Bp>Bp<Bp;Bp<Bp<Bp>Bp?Bp>Bp>Bp>Bp<Bp?Bp?Bp>Bp>Bp?Bp<Bp>Bp@BqDBqCBqABqABqBBqBBrGBrGBqBBqCBrEBrGBrGBrGBrHBrEBrIBrKBrKBrKBsPBsPBwkB{�B�B��B��B��B�B�B�2B�3B�3B�3B�:B�:B�AB�@B�NB�[B�qB�~B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�~B�~B�qB�qB�jB�hB�XB�TB�SB�SG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�B~�B��B~�Bu_BjB`�Bc�BjBl$Bn1Bp>BtWByuBz{B{�B|�B{�B{�B|�B~�B��B�B~�B}�B|�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�
B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�&B�(B�)B�)B�$B�#B�#B�$B�&B�$B�$B�#B�#B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B�B� B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��G�O�B��B��B��B��B��B|�BwoBxxB{�Bz�B��B�B�gB��B��B��B�B�B�B�B�411111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111144444444111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  G�O�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt( sw_cndr(PSAL,TEMP,PRES), TEMP, PRES_ADJUSTED )                                                                                                                                                                                         dP =-0.73 dbar.                                                                                                                                                                                                                                                 none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressures adjusted by using pressure offset at the sea surface. The quoted error is manufacturer specified accuracy in dbar.                                                                                                                                    The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   No significant salinity drift detected. The quoted error is max[0.01, 1xOWC uncertainty] in PSS-78.                                                                                                                                                             201904051344012019040513440120190405134401  AO  ARCAADJP                                                                    20181006004028    IP                G�O�G�O�G�O�                AO  ARGQQCPL                                                                    20181006004028  QCP$                G�O�G�O�G�O�5F03E           AO  ARGQQCPL                                                                    20181006004028  QCF$                G�O�G�O�G�O�C000            UW  ARSQUWQC    WOD & nearby Argo as visual check                               20190405134401  IP                  G�O�G�O�G�O�                