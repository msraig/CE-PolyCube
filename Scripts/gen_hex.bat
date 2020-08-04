@echo off
set DATA_ROOT_PATH=D:\data\testhexdata
set HEXEX_PATH=D:\Download\HexEx_Windows_1_01\HexEx_Windows
set HEX_OPT_PATH=D:\code\code_fromLiu\HexOpt\StandAlone\Frame2Mesh\Bin\msvc64
set POLYCUT_PATH=D:\polycut_public_release_20180929

set WORKING_PATH=%DATA_ROOT_PATH%\%1
echo DEFORMATION
..\Bin\SplitTet.exe %WORKING_PATH%\%1.vtk %WORKING_PATH%\%1_split.vtk > %1_log.txt
..\Bin\vtk2obj.exe %WORKING_PATH%\%1_split.vtk %WORKING_PATH%\%1.obj >> %1_log.txt
..\Bin\vtk2ovm.exe %WORKING_PATH%\%1_split.vtk %WORKING_PATH%\%1_split.ovm >> %1_log.txt
..\Bin\ScissorPoly.exe -input_name %WORKING_PATH%\%1_split.ovm -input_fea_name %WORKING_PATH%\%1.fea -output_name %WORKING_PATH%\%1_deform.vtk -mode 0 >> %1_log.txt

echo LABELING
..\Bin\tetvtk2mesh.exe %WORKING_PATH%\%1_deform.vtk %WORKING_PATH%\%1_deform.mesh
%POLYCUT_PATH%\mesh2vtu.exe %WORKING_PATH%\%1_deform.mesh %WORKING_PATH%\%1_deform.vtu
%POLYCUT_PATH%\polycut.exe        %WORKING_PATH%\%1_deform.vtu                      %WORKING_PATH%\%1_segm.vtu                        3 >> %1_log.txt 2>&1
%POLYCUT_PATH%\cusy2.exe          %WORKING_PATH%\%1_segm.vtu                 %WORKING_PATH%\%1_segm_deform.vtu                   >> %1_log.txt 2>&1
%POLYCUT_PATH%\optimizer.exe      %WORKING_PATH%\%1_segm_deform.vtu          %WORKING_PATH%\%1_segm_deform_opt.vtu           100 >> %1_log.txt 2>&1
..\Bin\vtudecode.exe %WORKING_PATH%\%1_segm_deform_opt.vtu %WORKING_PATH%\%1_segm_deform_opt_decode.vtu >> %1_log.txt
..\Bin\vtu2vtklabel.exe %WORKING_PATH%\%1_segm_deform_opt_decode.vtu %WORKING_PATH%\%1_polycube.vtk %WORKING_PATH%\%1_chartlabel.txt >> %1_log.txt
del %WORKING_PATH%\%1_segm_deform_opt_decode.vtu %WORKING_PATH%\%1_polycube.vtk >> %1_log.txt

echo FLATTENING
..\Bin\ScissorPoly.exe -input_name %WORKING_PATH%\%1_deform.vtk -input_cl_name %WORKING_PATH%\%1_chartlabel.txt -output_name %WORKING_PATH%\%1_deform_polycube.vtk -mode 1 -sr 2 -smd 0 -ac 0 -hex_flag 1 >> %1_log.txt

echo FEATURE_TRANSFER
..\Bin\tetovm2vtk.exe %WORKING_PATH%\%1_split.ovm %WORKING_PATH%\%1_split.vtk >> %1_log.txt
..\Bin\ScissorPoly.exe -mode 2 -rcl 0 -input_name %WORKING_PATH%\%1_split.vtk -input_pq_name %WORKING_PATH%\%1_deform_polycube.vtk -input_cl_name %WORKING_PATH%\%1_chartlabel.txt -input_fea_name %WORKING_PATH%\%1.fea -output_name %WORKING_PATH%\%1_polycube.tfe >> %1_log.txt

echo GENERATE_CUT
..\Bin\ScissorPoly.exe -mode 3 -bs 1 -cl 0 -mv 2.0 -ct 1 -sr 2 -th 1.0 -mr 0.5 -ac 1 -hex_flag 1 -fos 0.000001 -dcf 0 -input_name %WORKING_PATH%\%1_split.vtk -input_pq_name %WORKING_PATH%\%1_deform_polycube.vtk -input_cl_name %WORKING_PATH%\%1_chartlabel.txt -input_fea_name %WORKING_PATH%\%1_polycube.fea -output_name %WORKING_PATH%\%1_cut.vtk >> %1_log.txt

echo GENERATE_HEX
%HEXEX_PATH%\hexex.exe %WORKING_PATH%\%1_cut_ori.hexex %WORKING_PATH%\%1_cut_ori_hex.ovm >> %1_log.txt
%HEX_OPT_PATH%\HexaConvert.exe "-i" %WORKING_PATH%\%1_cut_ori_hex.ovm "-o" %WORKING_PATH%\%1_cut_ori_hex.vtk >> %1_log.txt
%HEXEX_PATH%\hexex.exe %WORKING_PATH%\%1_cut_polycube.hexex %WORKING_PATH%\%1_cut_polycube_hex.ovm >> %1_log.txt
%HEX_OPT_PATH%\HexaConvert.exe "-i" %WORKING_PATH%\%1_cut_polycube_hex.ovm "-o" %WORKING_PATH%\%1_cut_polycube_hex.vtk >> %1_log.txt

echo OPTIMIZE_HEX
..\Bin\TetHexFeatureTransfer.exe %WORKING_PATH%\%1_cut_polycube_hex.vtk %WORKING_PATH%\%1_polycube_segm.psfe %WORKING_PATH%\%1_cut_ori_hex.vtk %WORKING_PATH%\%1_cut_hex.hfe >> %1_log.txt
%HEX_OPT_PATH%\HexaConvert.exe -i %WORKING_PATH%\%1_cut_ori_hex.vtk -o %WORKING_PATH%\%1_cut_ori_hex_merge.vtk -m -p 1.0e-7 >> %1_log.txt
%HEX_OPT_PATH%\HexaRefiner.exe -e -t 0.4 -i %WORKING_PATH%\%1_cut_ori_hex_merge.vtk -o %WORKING_PATH%\%1_hex_opt.vtk -b %WORKING_PATH%\%1.obj -k %WORKING_PATH%\%1_polycube.tfe -l %WORKING_PATH%\%1_cut_hex.hfe >> %1_log.txt
del %WORKING_PATH%\%1_cut_ori_hex_merge.vtk %WORKING_PATH%\%1_cut_hex.hfe %WORKING_PATH%\%1_cut_polycube_hex.vtk %WORKING_PATH%\%1_cut_polycube_hex.ovm %WORKING_PATH%\%1_cut_ori_hex.vtk %WORKING_PATH%\%1_cut_ori_hex.ovm %WORKING_PATH%\%1_cut_polycube.hexex %WORKING_PATH%\%1_cut_ori.hexex
del %WORKING_PATH%\%1_polycube_segm.psfe %WORKING_PATH%\%1_polycube.fea %WORKING_PATH%\%1_polycube.tfe %WORKING_PATH%\%1_cut_ori.vtk %WORKING_PATH%\%1_split.vtk %WORKING_PATH%\%1_chartlabel.txt %WORKING_PATH%\%1_segm_deform_opt.vtu %WORKING_PATH%\segmentation_XYZ.obj %WORKING_PATH%\segmentation.mtl %WORKING_PATH%\%1_deform.vtu %WORKING_PATH%\%1_deform.mesh %WORKING_PATH%\%1_deform.fea %WORKING_PATH%\%1_deform.vtk %WORKING_PATH%\%1_split.ovm %WORKING_PATH%\%1_split.vtk %WORKING_PATH%\%1.obj
