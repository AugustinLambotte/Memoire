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
_FillValue                    �PArgo profile    3.1 1.2 19500101000000  20181006004028  20190405134401  5904772 US ARGO PROJECT                                                 STEPHEN RISER                                                   PRES            TEMP            PSAL               @A   AO  6627                            2C  D   APEX                            7392                            062915                          846 @�J�0*�w1   @�J��$&|@O&�t��A#t�j~�1   GPS     Primary sampling: mixed [deeper than nominal 985dbar: discrete; nominal 985dbar to surface: 2dbar-bin averaged]                                                                                                                                                    @A   B   B   @���@���@���A   A@  A`  A�  A�  A�  A�  A�  A�  A�  A�  A�33B  B  B  B   B(  B0  B8  B@  BH  BP  BX  B_��Bh  Bp  Bx  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B���B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  B�  C   C  C  C  C  C
  C  C  C  C  C  C  C  C  C  C  C   C"  C$  C&  C(  C*  C,  C.  C0  C2  C4  C6  C8  C:  C<  C>  C@  CB  CD  CF  CH  CJ  CL  CN  CP  CR  CT  CV  CX  CY�fC\  C^  C`  Cb  Cd  Cf  Ch  Cj  Cl  Cn  Cp  Cr  Ct  Cv  Cx  Cz  C|  C~  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��3C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C��C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  C�  D   D � D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D	  D	� D
  D
� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D  D� D   D � D!  D!� D"  D"� D#  D#� D$  D$� D%  D%� D&  D&� D'  D'� D(  D(� D)  D)� D*  D*� D+  D+� D,  D,� D-  D-� D.  D.� D/fD/� D0  D0� D1  D1� D2  D2� D3  D3� D4  D4� D5  D5� D6  D6� D7  D7� D8  D8� D9  D9� D:  D:� D;  D;� D<  D<� D=  D=� D>  D>� D?  D?� D@  D@� DA  DA� DB  DB� DC  DC� DD  DD� DE  DE� DF  DF� DG  DG� DH  DH� DI  DI� DJ  DJ� DK  DK� DL  DL� DM  DM� DN  DN� DO  DO� DP  DP� DQ  DQ� DR  DR� DS  DS� DT  DT� DU  DU� DV  DVy�DW  DW� DX  DX� DY  DY�fDZ  DZ� D[  D[� D\  D\� D]  D]� D^  D^� D_  D_� D`  D`� Da  Da� Db  Db� Dc  Dc� Dd  Dd� De  De� Df  Df� Dg  Dg� Dh  Dh� Di  Di� Dj  Dj� Dk  Dk� Dl  Dl� Dm  Dm� Dn  Dn� Do  Do� Dp  Dp� Dq  Dq� Dr  Dr� Ds  Ds� Dt  Dt� Dt�fDy�fD� D�33D���D��fD�3D�I�D�|�D�ٚD�3D�)�D�� DǼ�D��3D�S3Dڙ�D�� D��D�6fD� D�Ff111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @���@�(�A
{A+�AK�Ak�A��
A��
A��
A��
A��
A��
A��
A��
B�B
�B�B�B"�B*�B2�B:�BB�BJ�BR�BZ�Bb�Bj�Br�Bz�B�u�B�u�B�B�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�B�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�B�u�C ��C��C��C��C��C
��C��C��C��C��C��C��C��C��C��C��C ��C"��C$��C&��C(��C*��C,��C.��C0��C2��C4��C6��C8��C:��C<��C>��C@��CB��CD��CF��CH��CJ��CL��CN��CP��CR��CT��CV��CX��CZ�GC\��C^��C`��Cb��Cd��Cf��Ch��Cj��Cl��Cn��Cp��Cr��Ct��Cv��Cx��Cz��C|��C~��C�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�P�C�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�j>C�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qC�]qD .�D ��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D	.�D	��D
.�D
��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D.�D��D .�D ��D!.�D!��D".�D"��D#.�D#��D$.�D$��D%.�D%��D&.�D&��D'.�D'��D(.�D(��D).�D)��D*.�D*��D+.�D+��D,.�D,��D-.�D-��D..�D.��D/5D/��D0.�D0��D1.�D1��D2.�D2��D3.�D3��D4.�D4��D5.�D5��D6.�D6��D7.�D7��D8.�D8��D9.�D9��D:.�D:��D;.�D;��D<.�D<��D=.�D=��D>.�D>��D?.�D?��D@.�D@��DA.�DA��DB.�DB��DC.�DC��DD.�DD��DE.�DE��DF.�DF��DG.�DG��DH.�DH��DI.�DI��DJ.�DJ��DK.�DK��DL.�DL��DM.�DM��DN.�DN��DO.�DO��DP.�DP��DQ.�DQ��DR.�DR��DS.�DS��DT.�DT��DU.�DU��DV.�DV�RDW.�DW��DX.�DX��DY.�DY�DZ.�DZ��D[.�D[��D\.�D\��D].�D]��D^.�D^��D_.�D_��D`.�D`��Da.�Da��Db.�Db��Dc.�Dc��Dd.�Dd��De.�De��Df.�Df��Dg.�Dg��Dh.�Dh��Di.�Di��Dj.�Dj��Dk.�Dk��Dl.�Dl��Dm.�Dm��Dn.�Dn��Do.�Do��Dp.�Dp��Dq.�Dq��Dr.�Dr��Ds.�Ds��Dt.�Dt��Dt�Dy�D�'\D�J�D��)D���D��D�`�D��)D���D��D�@�D��\D��)D�
�D�j�Dڰ�D��\D�$)D�M�D�\D�]�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111  @��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@�E�@�M�@�M�@�V@���@��!@���@�~�@���@�
=@��y@�ȴ@���@�;d@�33@�33@�33@�33@�33@�33@�33@�;d@�;d@�C�@�;d@�33@�33@�S�@�S�@�K�@�C�@�33@�33@�33@�;d@�C�@�S�@�K�@�;d@�;d@�\)@�C�@�C�@�K�@�K�@�K�@�S�@�K�@�K�@�\)@�|�@�|�@�|�@��@���@��m@���@��P@��F@�b@��@�|�@���@� �@� �@� �@� �@� �@� �@��@��@�b@�b@��@��@��@�b@� �@� �@�(�@�1'@�1'@�1'@�1'@�9X@�9X@�9X@� �@��@� �@�(�@�1'@�1@�1'@�I�@�9X@�1'@�(�@� �@�(�@��@��m@�  @�  @��m@��;@���@�b@�b@��@� �@�b@�1@�  @�1@�1@���@���@���@���@��@���@�1@���@��@��@��@��m@��m@���@�  @�1@�1@�1@�b@�b@�b@�b@�b@��@��@��@�b@��@�b@� �@�  @�b@� �@�A�@��w@��P@�K�@��H@���@���@��\@��\@�n�@�V@�-@��@�@�p�@�G�@���@��D@�bN@�Ĝ@�n�@��y@�Q�@��@��
@�1@�Z@��D@��D@���@���@�Ĝ@��D@��D@��@�j@�A�@��@�ƨ@��@�dZ@�33@�@��y@���@��!@���@��\@�n�@�^5@�^5@�E�@�-@��@�p�@�G�@�7L@�V@��@�Ĝ@���@��D@��D@��D@�j@�A�@��@�w@��@~�@~@}�@}�@}�T@}�h@|��@}V@|1@{��@z�@zJ@y��@yX@x��@xr�@xbN@x1'@w��@w�w@wl�@v��@v��@v�R@v�R@v��@vv�@vV@vV@vV@vV@v$�@v{@v@u�-@u�@up�@u�@up�@u�@t�j@t�D@tI�@s��@s��@s��@t1@t1@s��@s�F@s�
@s�m@s�F@sƨ@sƨ@s��@s��@sdZ@s�@s�@s33@s33@s"�@s"�@r�H@rn�@rM�@r=q@rJ@q�^@q�7@q7L@qG�@q�@q&�@q�@qG�@q�@q%@p�u@p�u@pQ�@p1'@pA�@pA�@p �@o�@o��@o��@o�P@ol�@o�@n�@n$�@m�@m�T@m�-@mp�@m?}@m/@m/@m?}@m?}@m/@l��@l�/@l��@l�j@l�@l��@l��@l��@l�D@l�D@lz�@lj@lZ@lZ@lZ@lI�@lI�@l9X@l�@l1@l1@l1@l1@l1@k��@k��@k�m@k�m@kƨ@k��@k�@kS�@k33@ko@k@j�H@j��@j��@j��@j�\@j~�@jn�@j^5@j-@j�@i��@i�@i�#@i�^@i�^@i�^@i�^@i�^@i��@i��@i��@ix�@ihs@iX@iG�@i7L@i�@i�@h��@h�`@hĜ@h�9@h��@h�@h�@hr�@hQ�@hQ�@hA�@hA�@hA�@h1'@h �@hb@h  @h  @g�;@g�;@g�;@g�w@g��@g��@g�P@g�P@g��@g�P@gl�@gK�@g+@g�@g
=@f��@f�y@f�@f�R@f��@f��@f��@f�+@fv�@fff@fff@fE�@f{@e�@e�@e��@e@e�-@e�-@e�-@e��@e�h@e�@e`B@eO�@e?}@e/@e?}@e?}@e/@e�@d��@d�@d�/@d�j@d�D@dz�@dj@dI�@d9X@dI�@dZ@dI�@d9X@d�@c��@c�m@cƨ@c��@c�m@d1@c�m@cƨ@c��@cdZ@cS�@cS�@cS�@cS�@c33@c@b�@b�H@b�H@b�H@b�!@b^5@b^5@bn�@bn�@bJ@a�7@_�;@^�R@]p�@\�/@[��@[@Z��@^��@a��@e��@hĜ@ix�@ix�@h1'@f��@d��@aG�@]�-@Y�^@VV111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  @�E�@�M�@�M�@�V@���@��!@���@�~�@���@�
=@��y@�ȴ@���@�;d@�33@�33@�33@�33@�33@�33@�33@�;d@�;d@�C�@�;d@�33@�33@�S�@�S�@�K�@�C�@�33@�33@�33@�;d@�C�@�S�@�K�@�;d@�;d@�\)@�C�@�C�@�K�@�K�@�K�@�S�@�K�@�K�@�\)@�|�@�|�@�|�@��@���@��m@���@��P@��F@�b@��@�|�@���@� �@� �@� �@� �@� �@� �@��@��@�b@�b@��@��@��@�b@� �@� �@�(�@�1'@�1'@�1'@�1'@�9X@�9X@�9X@� �@��@� �@�(�@�1'@�1@�1'@�I�@�9X@�1'@�(�@� �@�(�@��@��m@�  @�  @��m@��;@���@�b@�b@��@� �@�b@�1@�  @�1@�1@���@���@���@���@��@���@�1@���@��@��@��@��m@��m@���@�  @�1@�1@�1@�b@�b@�b@�b@�b@��@��@��@�b@��@�b@� �@�  @�b@� �@�A�@��w@��P@�K�@��H@���@���@��\@��\@�n�@�V@�-@��@�@�p�@�G�@���@��D@�bN@�Ĝ@�n�@��y@�Q�@��@��
@�1@�Z@��D@��D@���@���@�Ĝ@��D@��D@��@�j@�A�@��@�ƨ@��@�dZ@�33@�@��y@���@��!@���@��\@�n�@�^5@�^5@�E�@�-@��@�p�@�G�@�7L@�V@��@�Ĝ@���@��D@��D@��D@�j@�A�@��@�w@��@~�@~@}�@}�@}�T@}�h@|��@}V@|1@{��@z�@zJ@y��@yX@x��@xr�@xbN@x1'@w��@w�w@wl�@v��@v��@v�R@v�R@v��@vv�@vV@vV@vV@vV@v$�@v{@v@u�-@u�@up�@u�@up�@u�@t�j@t�D@tI�@s��@s��@s��@t1@t1@s��@s�F@s�
@s�m@s�F@sƨ@sƨ@s��@s��@sdZ@s�@s�@s33@s33@s"�@s"�@r�H@rn�@rM�@r=q@rJ@q�^@q�7@q7L@qG�@q�@q&�@q�@qG�@q�@q%@p�u@p�u@pQ�@p1'@pA�@pA�@p �@o�@o��@o��@o�P@ol�@o�@n�@n$�@m�@m�T@m�-@mp�@m?}@m/@m/@m?}@m?}@m/@l��@l�/@l��@l�j@l�@l��@l��@l��@l�D@l�D@lz�@lj@lZ@lZ@lZ@lI�@lI�@l9X@l�@l1@l1@l1@l1@l1@k��@k��@k�m@k�m@kƨ@k��@k�@kS�@k33@ko@k@j�H@j��@j��@j��@j�\@j~�@jn�@j^5@j-@j�@i��@i�@i�#@i�^@i�^@i�^@i�^@i�^@i��@i��@i��@ix�@ihs@iX@iG�@i7L@i�@i�@h��@h�`@hĜ@h�9@h��@h�@h�@hr�@hQ�@hQ�@hA�@hA�@hA�@h1'@h �@hb@h  @h  @g�;@g�;@g�;@g�w@g��@g��@g�P@g�P@g��@g�P@gl�@gK�@g+@g�@g
=@f��@f�y@f�@f�R@f��@f��@f��@f�+@fv�@fff@fff@fE�@f{@e�@e�@e��@e@e�-@e�-@e�-@e��@e�h@e�@e`B@eO�@e?}@e/@e?}@e?}@e/@e�@d��@d�@d�/@d�j@d�D@dz�@dj@dI�@d9X@dI�@dZ@dI�@d9X@d�@c��@c�m@cƨ@c��@c�m@d1@c�m@cƨ@c��@cdZ@cS�@cS�@cS�@cS�@c33@c@b�@b�H@b�H@b�H@b�!@b^5@b^5@bn�@bn�G�O�@a�7@_�;@^�R@]p�@\�/@[��@[@Z��@^��@a��@e��@hĜ@ix�@ix�@h1'@f��@d��@aG�@]�-@Y�^@VV111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  ;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB� B� B� B� B� B� B� B�B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B�B�B� B� B�B�B� B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B� B�B�B� B�B�B�B�B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B� B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�%B�%B�%B�%B�%B�+B�+B�7B�7B�7B�DB�DB�DB�PB�{B��B��B��B��B��B��B��B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�{B�{B�{B�{B�{B�{B�{B�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�uB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�oB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�hB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�bB�\B�\B�\B�\B�\B�\B�\B�\B�\B�\B�VB�VB�VB�VB�VB�VB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�PB�JB�JB�JB�JB�JB�JB�JB�JB�JB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�DB�=B�=B�=B�=B�DB�DB�=B�=B�=B�=B�=B�=B�7B�7B�7B�7B�7B�7B�7B�1B�1B�1B�1B�1B�JB�%B�B�B~�B~�B{�B|�B|�B�7B�uB��B�B�?B�RB�jB�qB�wB�jB�jB�dB�j111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  B�B�B�B�B�B�B�B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B�B�B��B��B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B��B��B�B��B��B��B��B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�B�,B�QB�QB�PB�^B�pB�sB�xB�}B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B�~B�~B��B�~B�B�|B�|B�zB�wB�vB�vB�rB�tB�mB�eB�gB�eB�eB�eB�`B�`B�`B�`B�XB�XB�WB�XB�XB�RB�QB�PB�WB�XB�XB�PB�YB�RB�RB�RB�SB�PB�RB�RB�MB�KB�KB�MB�MB�MB�MB�MB�KB�KB�MB�KB�KB�KB�MB�MB�KB�MB�KB�KB�KB�KB�EB�EB�FB�IB�FB�GB�HB�?B�BB�DB�@B�@B�@B�@B�BB�BB�CB�<B�<B�<B�>B�<B�<B�<B�4B�6B�4B�4B�4B�4B�-B�-B�(B�)B�'B�)B�#B�#B�#B�!B�#B�#B�!B�&B�#B�!B�#B�$B�#B�$B�!B�!B�#B�$B�$B�#B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�	B�	B�	B�B�	B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B�B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��B��G�O�B��B��B��B~�B~�B{�B|�B|�B��B�B�mB��B��B��B�B�B�B�B�B�
B�111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111114111111111111111111111  <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
PRES            TEMP            PSAL            PRES_ADJUSTED = PRES - dP                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            PSAL_ADJUSTED = sw_salt( sw_cndr(PSAL,TEMP,PRES), TEMP, PRES_ADJUSTED )                                                                                                                                                                                         dP =-0.73 dbar.                                                                                                                                                                                                                                                 none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Pressures adjusted by using pressure offset at the sea surface. The quoted error is manufacturer specified accuracy in dbar.                                                                                                                                    The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.                                                                                                                                                   No significant salinity drift detected. The quoted error is max[0.01, 1xOWC uncertainty] in PSS-78.                                                                                                                                                             201904051344012019040513440120190405134401  AO  ARCAADJP                                                                    20181006004028    IP                G�O�G�O�G�O�                AO  ARGQQCPL                                                                    20181006004028  QCP$                G�O�G�O�G�O�5F03E           AO  ARGQQCPL                                                                    20181006004028  QCF$                G�O�G�O�G�O�8000            UW  ARSQUWQC    WOD & nearby Argo as visual check                               20190405134401  IP                  G�O�G�O�G�O�                