CDF      
      	DATE_TIME         	STRING256         STRING64   @   STRING32       STRING16      STRING8       STRING4       STRING2       N_PROF        N_PARAM       N_LEVELS   \   	N_HISTORY          N_CALIB          
   title         Argo float vertical profile    institution       CORIOLIS   source        
Argo float     history       S2022-08-04T07:30:39Z creation; 2022-08-25T13:29:24Z last update (BSH ARSQ software)    
references        (http://www.argodatamgt.org/Documentation   user_manual_version       3.1    Conventions       Argo-3.1 CF-1.6    featureType       trajectoryProfile      decoder_version       	CODA_050n      comment_dmqc_operator         CPRIMARY | https://orcid.org/0000-0003-2129-3325 | Birgit Klein, BSH       C   	DATA_TYPE                  	long_name         	Data type      conventions       Argo reference table 1     
_FillValue                    9�   FORMAT_VERSION                 	long_name         File format version    
_FillValue                    9�   HANDBOOK_VERSION               	long_name         Data handbook version      
_FillValue                    :    REFERENCE_DATE_TIME                 	long_name         !Date of reference for Julian days      conventions       YYYYMMDDHHMISS     
_FillValue                    :   DATE_CREATION                   	long_name         Date of file creation      conventions       YYYYMMDDHHMISS     
_FillValue                    :   DATE_UPDATE                 	long_name         Date of update of this file    conventions       YYYYMMDDHHMISS     
_FillValue                    :$   PLATFORM_NUMBER                   	long_name         Float unique identifier    conventions       WMO float identifier : A9IIIII     
_FillValue                    :4   PROJECT_NAME                  	long_name         Name of the project    
_FillValue                  @  :<   PI_NAME                   	long_name         "Name of the principal investigator     
_FillValue                  @  :|   STATION_PARAMETERS           	            	long_name         ,List of available parameters for the station   conventions       Argo reference table 3     
_FillValue                  @  :�   CYCLE_NUMBER               	long_name         Float cycle number     conventions       =0...N, 0 : launch cycle (if exists), 1 : first complete cycle      
_FillValue         ��        :�   	DIRECTION                  	long_name         !Direction of the station profiles      conventions       -A: ascending profiles, D: descending profiles      
_FillValue                    ;    DATA_CENTRE                   	long_name         .Data centre in charge of float data processing     conventions       Argo reference table 4     
_FillValue                    ;   DC_REFERENCE                  	long_name         (Station unique identifier in data centre   conventions       Data centre convention     
_FillValue                     ;   DATA_STATE_INDICATOR                  	long_name         1Degree of processing the data have passed through      conventions       Argo reference table 6     
_FillValue                    ;(   	DATA_MODE                  	long_name         Delayed mode or real time data     conventions       >R : real time; D : delayed mode; A : real time with adjustment     
_FillValue                    ;,   PLATFORM_TYPE                     	long_name         Type of float      conventions       Argo reference table 23    
_FillValue                     ;0   FLOAT_SERIAL_NO                   	long_name         Serial number of the float     
_FillValue                     ;P   FIRMWARE_VERSION                  	long_name         Instrument firmware version    
_FillValue                     ;p   WMO_INST_TYPE                     	long_name         Coded instrument type      conventions       Argo reference table 8     
_FillValue                    ;�   JULD               	long_name         ?Julian day (UTC) of the station relative to REFERENCE_DATE_TIME    standard_name         time   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~       axis      T           ;�   JULD_QC                	long_name         Quality on date and time   conventions       Argo reference table 2     
_FillValue                    ;�   JULD_LOCATION                  	long_name         @Julian day (UTC) of the location relative to REFERENCE_DATE_TIME   units         "days since 1950-01-01 00:00:00 UTC     conventions       8Relative julian days with decimal part (as parts of day)   
resolution        >�EȠ�Q)   
_FillValue        A.�~            ;�   LATITUDE               	long_name         &Latitude of the station, best estimate     standard_name         latitude   units         degree_north   
_FillValue        @�i�       	valid_min         �V�        	valid_max         @V�        axis      Y           ;�   	LONGITUDE                  	long_name         'Longitude of the station, best estimate    standard_name         	longitude      units         degree_east    
_FillValue        @�i�       	valid_min         �f�        	valid_max         @f�        axis      X           ;�   POSITION_QC                	long_name         ,Quality on position (latitude and longitude)   conventions       Argo reference table 2     
_FillValue                    ;�   POSITIONING_SYSTEM                    	long_name         Positioning system     
_FillValue                    ;�   PROFILE_MTIME_QC               	long_name         $Global quality flag of MTIME profile   conventions       Argo reference table 2a    
_FillValue                    ;�   PROFILE_PRES_QC                	long_name         #Global quality flag of PRES profile    conventions       Argo reference table 2a    
_FillValue                    ;�   PROFILE_TEMP_QC                	long_name         #Global quality flag of TEMP profile    conventions       Argo reference table 2a    
_FillValue                    ;�   PROFILE_PSAL_QC                	long_name         #Global quality flag of PSAL profile    conventions       Argo reference table 2a    
_FillValue                    ;�   VERTICAL_SAMPLING_SCHEME                  	long_name         Vertical sampling scheme   conventions       Argo reference table 16    
_FillValue                    ;�   CONFIG_MISSION_NUMBER                  	long_name         :Unique number denoting the missions performed by the float     conventions       !1...N, 1 : first complete mission      
_FillValue         ��        <�   MTIME            
         	long_name         LFractional day of the individual measurement relative to JULD of the station   
_FillValue        A.�~       units         days   	valid_min         �         	valid_max         @         C_format      %.6f   FORTRAN_format        F.6    
resolution        5�7�     �  <�   MTIME_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  \  ?�   PRES         
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        p  @   PRES_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  \  A�   PRES_ADJUSTED            
      
   	long_name         )Sea water pressure, equals 0 at sea-level      standard_name         sea_water_pressure     
_FillValue        G�O�   units         decibar    	valid_min                	valid_max         F;�    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���   axis      Z        p  A�   PRES_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  \  CP   PRES_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         decibar    C_format      %7.1f      FORTRAN_format        F7.1   
resolution        =���     p  C�   TEMP         
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  E   TEMP_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  \  F�   TEMP_ADJUSTED            
      	   	long_name         $Sea temperature in-situ ITS-90 scale   standard_name         sea_water_temperature      
_FillValue        G�O�   units         degree_Celsius     	valid_min         �      	valid_max         B      C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  F�   TEMP_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  \  HX   TEMP_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         degree_Celsius     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  H�   PSAL         
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  J$   PSAL_QC          
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  \  K�   PSAL_ADJUSTED            
      	   	long_name         Practical salinity     standard_name         sea_water_salinity     
_FillValue        G�O�   units         psu    	valid_min         @      	valid_max         B$     C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  K�   PSAL_ADJUSTED_QC         
         	long_name         quality flag   conventions       Argo reference table 2     
_FillValue                  \  M`   PSAL_ADJUSTED_ERROR          
         	long_name         VContains the error on the adjusted values as determined by the delayed mode QC process     
_FillValue        G�O�   units         psu    C_format      %9.3f      FORTRAN_format        F9.3   
resolution        :�o     p  M�   HISTORY_INSTITUTION                      	long_name         "Institution which performed action     conventions       Argo reference table 4     
_FillValue                    [�   HISTORY_STEP                     	long_name         Step in data processing    conventions       Argo reference table 12    
_FillValue                    [�   HISTORY_SOFTWARE                     	long_name         'Name of software which performed action    conventions       Institution dependent      
_FillValue                    [�   HISTORY_SOFTWARE_RELEASE                     	long_name         2Version/release of software which performed action     conventions       Institution dependent      
_FillValue                    [�   HISTORY_REFERENCE                        	long_name         Reference of database      conventions       Institution dependent      
_FillValue                  @  [�   HISTORY_DATE                      	long_name         #Date the history record was created    conventions       YYYYMMDDHHMISS     
_FillValue                    [�   HISTORY_ACTION                       	long_name         Action performed on data   conventions       Argo reference table 7     
_FillValue                    \   HISTORY_PARAMETER                        	long_name         (Station parameter action is performed on   conventions       Argo reference table 3     
_FillValue                    \   HISTORY_START_PRES                    	long_name          Start pressure action applied on   
_FillValue        G�O�   units         decibar         \   HISTORY_STOP_PRES                     	long_name         Stop pressure action applied on    
_FillValue        G�O�   units         decibar         \   HISTORY_PREVIOUS_VALUE                    	long_name         +Parameter/Flag previous value before action    
_FillValue        G�O�        \    HISTORY_QCTEST                       	long_name         <Documentation of tests performed, tests failed (in hex form)   conventions       EWrite tests performed when ACTION=QCP$; tests failed when ACTION=QCF$      
_FillValue                    \$   	PARAMETER               	            	long_name         /List of parameters with calibration information    conventions       Argo reference table 3     
_FillValue                  @  O,   SCIENTIFIC_CALIB_EQUATION               	            	long_name         'Calibration equation for this parameter    
_FillValue                    Ol   SCIENTIFIC_CALIB_COEFFICIENT            	            	long_name         *Calibration coefficients for this equation     
_FillValue                    Sl   SCIENTIFIC_CALIB_COMMENT            	            	long_name         .Comment applying to this parameter calibration     
_FillValue                    Wl   SCIENTIFIC_CALIB_DATE               	             	long_name         Date of calibration    conventions       YYYYMMDDHHMISS     
_FillValue                  8  [lArgo profile    3.1 1.2 19500101000000  20220804073039  20220825132924  6900883 AWI                                                             Gerd ROHARDT                                                    MTIME           PRES            TEMP            PSAL               )A   IF                                  2C  D   NEMO                            220                             23-May-2012                     860 @֨T��+1   @֨T�� @S��<@��%�1   GPS         A   B   F   Primary sampling: discrete []                                                                                                                                                                                                                                      ����@  ��#Eh�  ��d�
�  ���m;   ����Z   ���O�   ��)���  ��_1��  ��	,@  ���w`  ��韬   ��+=@  �� �/@  ��m�6   ���d�   ���ޡ   ���GN   ��!�   ��z�@  ��{�v�  ��n]M@  ��"��  ��EȡP  ����߀  ���kU  ��
��`  ���%+P  ��g��p  ���Y �  �����  ��ۗ�  ���A<   ��:�   ���c�@  ��N �@  ����6P  ��0���  ��u0�   ���>��  ��A;��  ��˩��  ��.�  ��M^o�  �� a`  �����  ��X�%�  ���Ӏ  ����ݸ  ��So��  ��i�x  ���u1(  ��lw��  �� ���  ��؎�8  ����A�  ��!�@  ����P@  ������  ��<M^�  ���GM8  ���[g  ��G��  ��� �p  ���/p  ��o�   ����  ����5�  ���ò�  ���`�  ���8�  ��x���  ��RL�`  ��2�0  �����  ����uP  ����Q�  ����.  ������  ��`T�\  ��&�7�  �����  ��6  ��F��F  ��r8  �����  �ɉJ��  ��Er�  ��     ��ò��  ��A;�8  ��I��  ����F�  00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000=���?   @S33@�33@�  A  A)��AS33At��A�  A���A�  A���A���Aə�A�  A�33B33B��B��B��B$��B+33B4��B;��BA��BY33Bm��B~ffB���B���B�ffB�ffB���B���B�  B���B���C  C
��C�C  C(��C3�C<��CGL�CP�fCZ� Ce33Co33Cy  C�Y�C�Y�C�� C��fC��3C�Y�C��fC�� C�Y�C��3C�Y�C��C�ffC�� Cǌ�C�ٚC�� C�  C�ffD�D	33Dy�D�fD�D",�D(� D.�fD;@ DGٚD`s3Dy��D�ffD��3D�p D��D�\�D��3D�i�D��3D�l�D��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111=���?   @S33@�33@�  A  A)��AS33At��A�  A���A�  A���A���Aə�A�  A�33B33B��B��B��B$��B+33B4��B;��BA��BY33Bm��B~ffB���B���B�ffB�ffB���B���B�  B���B���C  C
��C�C  C(��C3�C<��CGL�CP�fCZ� Ce33Co33Cy  C�Y�C�Y�C�� C��fC��3C�Y�C��fC�� C�Y�C��3C�Y�C��C�ffC�� Cǌ�C�ٚC�� C�  C�ffD�D	33Dy�D�fD�D",�D(� D.�fD;@ DGٚD`s3Dy��D�ffD��3D�p D��D�\�D��3D�i�D��3D�l�D��11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@��@ύP@�|�@�o@�ȴ@���@�=q@�`B@�ƨ@�n�@�7L@��
@�o@�V@�S�@�{@�(�@�{@�;d@�dZ@�o@�Z@���@�@�`B@�Z@~��@xQ�@r�\@l�@g�@d�@ax�@]�@ZJ@V{@R-@Fȴ@7
=@,�D@#dZ@#�@"J@�R@-@��@
M�@��?�5??�1'?��?���?�J?��T?�S�?h1'?XQ�?8Q�?�-?�!?{?�7>���>��>���>�t�>��>#�
=���<�1=P�`;��
������;d���&�y���˾��u��Z��Q���Ͼ�vɿ	���P�vɿ&$ݿ,�D�2n��9���?;d�CS��F�y�G�11444444441111111111111111111111111111144111111111111111111111111111111111111111111111111111@ύP@�|�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�@��
@�o@�V@�S�@�{@�(�@�{@�;d@�dZ@�o@�Z@���@�@�`B@�Z@~��@xQ�@r�\@l�@g�@d�@ax�@]�@ZJ@V{@R-@Fȴ@7
=@,�DG�O�G�O�@"J@�R@-@��@
M�@��?�5??�1'?��?���?�J?��T?�S�?h1'?XQ�?8Q�?�-?�!?{?�7>���>��>���>�t�>��>#�
=���<�1=P�`;��
������;d���&�y���˾��u��Z��Q���Ͼ�vɿ	���P�vɿ&$ݿ,�D�2n��9���?;d�CS��F�y�G�11444444441111111111111111111111111111144111111111111111111111111111111111111111111111111111;o;oG�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oG�O�G�O�;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;o;oB�3B�!B��Bl�BgmB-B
�HB
}�B
>wB
PB	�`B	��B	}�B	]/B	VB	k�B	}�B	�B
s�B
�3B  B9XBJ�BP�BO�BO�BQ�BVBk�B~�B{�Br�BjBn�BdZBYB?}B�B�B%B�5B��B��Bz�B�BjBaHBk�BjB}�B��B��B�hB��B�LB��B�-B�HB�yB�B��B�HB�B�/B�B�NB��B�`BoB  B�NB��BB+B�B��BBB��B	7B	7B1B	7B
=B
=BDBDBPBJBJBVBh33444444444444444333333333333333333333344333333333333333333333333333333333333333333333333333G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�44444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�G�O�MTIME           PRES            TEMP            PSAL            Not applicable                                                                                                                                                                                                                                                  PRES_ADJUSTED (cycle i) = PRES (cycle i) - Surface Pressure (cycle i+1)                                                                                                                                                                                         TEMP_ADJUSTED = TEMP                                                                                                                                                                                                                                            PSAL_ADJUSTED = PSAL (re-calculated by using PRES_ADJUSTED)                                                                                                                                                                                                     Not applicable                                                                                                                                                                                                                                                  Surface pressure = 0 dbar                                                                                                                                                                                                                                       none                                                                                                                                                                                                                                                            none                                                                                                                                                                                                                                                            Not applicable                                                                                                                                                                                                                                                  No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar                                                                                                                                                            No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90                                                                                                                                          Severe sensor malfunction, measurements are bad and not adjustable                                                                                                                                                                                              20220825132924202208251329242022082513292420220825132924IF  ARFMCODA050n                                                                20220804073039                      G�O�G�O�G�O�                IF  ARGQCOQC5.8                                                                 20220804073217  QCP$                G�O�G�O�G�O�000000000208F37EIF  ARGQCOQC5.8                                                                 20220804073217  QCF$                G�O�G�O�G�O�000000000200C000GE  ARSQOW  1.0 ARGO CTD ref. database: CTD_for_DMQC_2021V02 + ARGO climatology 20220825132924  IP  PSAL            =���D��G�O�                