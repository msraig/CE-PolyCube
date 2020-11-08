#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/filesystem.hpp>
#include <gflags/gflags.h>
#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include "TetStructure.h"
#include "Polycube_Algorithm.h"
#include "SmallVec.h"
DEFINE_int32(mode, 3, "processing mode"); //0: deform initial object(ovm format), 1: flattening deformed object(vtk format), 2: feature transfer, 3: cut
DEFINE_string(input_name, "./input/cup_input.vtk", "input file name");
DEFINE_string(output_name, "output.vtk", "output file name");
DEFINE_string(input_cl_name, "input.txt", "input chart label name");
DEFINE_string(input_pq_name, "input_polycube.vtk", "input polycube name");
DEFINE_string(input_fea_name, "nofeature.fea", "input feature name");
DEFINE_string(log_dir, "./", "log directory");
DEFINE_int32(hex_flag, 1, "Hex Meshing Flag"); //to produce float result, hex flag should be set to 0; to produce int result, set it as 1
DEFINE_double(ss, 1.0, "sigma_s");
DEFINE_double(sr, 3.0, "sigma_r"); //controling the target edge length of hex mesh
DEFINE_double(cd, 0.3, "cut depth");
DEFINE_int32(dt, 0, "distortion type"); //0: isometric, 1: conformal 2: volumetric
DEFINE_double(mv, 2.0, "max volume");
DEFINE_int32(cf, 1, "cut flag");
DEFINE_int32(si, 1, "save intermediate result");
DEFINE_int32(cl, 0, "set cube length");
DEFINE_int32(bs, 5, "cut batch size");
DEFINE_double(cr, 1.0, "cut edge ratio");
DEFINE_int32(ii, 10, "ivf iteration");
DEFINE_int32(sf, 0, "save flip information");
DEFINE_double(mr, 0.5, "min depth ratio");
DEFINE_double(ld, 0.5, "lambda for flattening");
DEFINE_int32(ac, 0, "add constraints for flattening");
DEFINE_int32(ct, 0, "cut type"); // 0: auto cut 1: auto_cut_max 2: auto_cut_min 3: prescribe_cut
DEFINE_int32(cc, 0, "cut all convex edge"); //all convex edge in ignored_edges
DEFINE_int32(sh, 0, "sample hex");
DEFINE_int32(cn, -1, "number of cut");
DEFINE_int32(ipt, 0, "IVF polycube type"); // 0: IVF standard; 1: optimize polycube;
DEFINE_int32(smd, 0, "Set min difference for preprocessing"); //work only when hex_flag is set to 0, if used, right part of all constraints is cube_length.
DEFINE_int32(is, 0, "Integer solver for flattening"); //work only when hex_flag is set to 1, and when chart_value_ori is empty
DEFINE_double(th, 1.8, "threshold for polycube edge distortion");
DEFINE_int32(rd, 0, "Randomize cut edge");
DEFINE_int32(dff, 1, "Deformation For Flattening");
DEFINE_double(smooth_factor, 0.1, "Smooth Factor");
DEFINE_double(fsr, 0.8, "Feature Select Ratio");
DEFINE_int32(rcl, 1, "Repair Chart Label"); //used for mode 2
DEFINE_double(fos, 0.0, "Flattening Offset");
DEFINE_int32(dcf, 0, "Double Cube Length Flag"); //double the cube_length when calculating chart value

int main(int argc, char** argv)
{
	gflags::ParseCommandLineFlags(&argc, &argv, true);
	if (FLAGS_mode == 0)
	{
		//feature-based deform, only support ovm file
		std::string input_model_name = FLAGS_input_name;
		std::string input_fea_name = FLAGS_input_fea_name;
		std::string output_model_name = FLAGS_output_name;
		std::string output_fea_name = output_model_name.substr(0, output_model_name.length() - 3) + "fea";
		double sigma_r = FLAGS_sr;
		Polycube_Algorithm pa;
		VolumeMesh mesh_;
		OpenVolumeMesh::IO::FileManager fileManager;
		fileManager.readFile(input_model_name, mesh_);
		std::cout << "Mesh Info:\n" << "nv: " << mesh_.n_vertices() << " nc: " << mesh_.n_cells() << std::endl;
		pa.SetMesh(&mesh_);
		pa.SetSigmaR(sigma_r);
		bool feature_exist = false;
		if (pa.load_feature_edges(input_fea_name.c_str()))
		{
			feature_exist = true;
			//deform polyline
			bool polyline_guidance = true;
			bool ns_flag = false;
			double pq_facet_ratio_threshold = 0.08;
			pa.repair_features(polyline_guidance, ns_flag, FLAGS_smooth_factor, pq_facet_ratio_threshold);
		}
		//feature-aware deformation
		double pq_facet_ratio_threshold = 0.08;
		int deform_iter = 5;
		pa.feature_aware_deformation(deform_iter, 1.0, pq_facet_ratio_threshold);
		pa.save_tet_vtk(output_model_name.c_str());
		if (feature_exist)
		{
			pa.save_feature_edges_ovm(output_fea_name.c_str());
		}
	}
	else if (FLAGS_mode == 1)
	{
		//flattening deformed object
		std::string input_model_name = FLAGS_input_name;
		std::string input_cl_name = FLAGS_input_cl_name;
		std::string output_model_name = FLAGS_output_name;
		
		int last_slash = (int)input_model_name.find_last_of("\\") + 1;
		std::string mainname = input_model_name.substr(last_slash, input_model_name.length() - last_slash - 4);
		Polycube_Algorithm pa;
		TetStructure<double> input_ori;
		input_ori.load_vtk(input_model_name.c_str());
		if (input_ori.tetra_vertices.size() == 0) return 1;
		pa.SetMesh(&input_ori);
		pa.SetSigmaR(FLAGS_sr);
		if (!pa.load_chart_label(input_cl_name.c_str())) return 1;
		//additional deformation
		if (FLAGS_dff)
		{
			pa.do_deformation_constrained_vtk(FLAGS_ss);
			pa.do_deformation_constrained_continued();
		}
		pa.pq_flattening->set_hex_meshing_flag(FLAGS_hex_flag);
		//flattening
		pa.do_flattening_constrained_vtk(FLAGS_dt, -1.0, FLAGS_ld, FLAGS_ac, FLAGS_ipt, FLAGS_smd, FLAGS_is, FLAGS_fos, FLAGS_dcf);
		pa.save_tet_vtk(FLAGS_output_name.c_str());
	}
	else if (FLAGS_mode == 2)
	{
		//feature transfer
		//intput: 3d shape, chart label, feature
		//output: feature(including featured tc, index only) feature_long_vtk feature_short_vtk
		Polycube_Algorithm pa;
		TetStructure<double> input_ori;
		TetStructure<double> input_polycube;
		int last_slash = (int)FLAGS_input_name.find_last_of("\\") + 1;
		std::string mainname = FLAGS_input_name.substr(last_slash, FLAGS_input_name.length() - last_slash - 4);
		input_ori.load_vtk(FLAGS_input_name.c_str());
		input_polycube.load_vtk(FLAGS_input_pq_name.c_str());
		//pa.SetMesh(&input_ori);
		pa.SetMesh(&input_polycube);
		pa.SetSigmaR(FLAGS_sr);
		int flip_count = 0;
		const std::vector<Tetrahedron<double>*> &tetras_polycube = input_polycube.tetras;
		int nc = (int)tetras_polycube.size();
		for (size_t i = 0; i < nc; i++)
		{
			if (tetras_polycube[i]->compute_tet_volume() < 0)
				flip_count++;
		}
		if (!pa.load_chart_label(FLAGS_input_cl_name.c_str()))
		{
			std::cout << "Loading chart label error!" << std::endl;
			return 1;
		}
		if (!pa.load_feature_edges(FLAGS_input_fea_name.c_str()))
		{
			std::cout << "Loading Feature error!" << std::endl;
			return 1;
		}
		pa.flattening_prepare_vtk();
		pa.pq_flattening->load_feature_edges_vtk(pa.init_feature_edge_array);
		int min_feature_length = 3;
		bool feature_extract_success = pa.pq_flattening->feature_extraction(FLAGS_fsr, min_feature_length); //first repair then repair would be better
		if (!feature_extract_success)
		{
			if (FLAGS_rcl)
			{
				//repair chart label
				bool corner_correct = pa.pq_flattening->repair_chartlabel_feature(8, 20);
				if (!corner_correct)
				{
					pa.pq_flattening->save_error_feature_ptid_ptsformat((FLAGS_output_name.substr(0, FLAGS_output_name.length() - 4) + ".pts").c_str());
				}
				pa.pq_flattening->save_chartlabel((FLAGS_output_name.substr(0, FLAGS_output_name.length() - 4) + "_chartlabel_repair.txt").c_str());
			}
		}
		pa.pq_flattening->save_feature_edges_vtk_tfeformat(FLAGS_output_name.c_str(), &input_ori);
		pa.pq_flattening->save_feature_edges_vtk_feaformat((FLAGS_output_name.substr(0, FLAGS_output_name.length() - 4) + ".fea").c_str());
		pa.pq_flattening->save_feature_long_edges_vtk_psfeformat((FLAGS_output_name.substr(0, FLAGS_output_name.length() - 4) + "_segm.psfe").c_str());
	}
	else if (FLAGS_mode == 3)
	{
		//cut new version
		//std::string output_time_prefix("./time/");
		//std::string output_flips_prefix("./flips/");
		//std::string output_intermediate_prefix("./intermediate/");
		std::string output_cubelength_prefix("./cubelength_cut_new/");
		//std::string output_cutnumber_prefix("./cutnum/");
		std::string input_cubelength_prefix("./cubelength/");
		//std::string output_flipinfo_prefix("./flipinfo/");
		std::string log_string;
		long start_t = clock();
		std::cout << "Scissor Polycube begin: " << std::endl;
		//created folder
		std::string input_ori_name = FLAGS_input_name;
		std::string output_name = FLAGS_output_name;
		std::string name_prefix = output_name.substr(0, output_name.length() - 4);
		std::string filename = output_name.substr(output_name.find_last_of('\\') + 1, output_name.find_last_of('.') - output_name.find_last_of('\\') - 1);
		std::string input_polycube_name = FLAGS_input_pq_name;
		//std::cout << "Input File2 : " << argv[2] << std::endl;
		std::string input_chart_label_name = FLAGS_input_cl_name;
		std::string input_feature_name = FLAGS_input_fea_name;
		std::cout << "Input origin mesh: " << input_ori_name << std::endl;
		std::cout << "Input polycube: " << input_polycube_name << std::endl;
		std::cout << "Input chartlabel: " << input_chart_label_name << std::endl;
		std::cout << "Input feature: " << input_feature_name << std::endl;
		std::cout << "file name: " << filename << std::endl;
		std::cout << "name prefix: " << name_prefix << std::endl;
		log_string += "Input origin mesh: " + input_ori_name + "\n";
		log_string += "Input polycube: " + input_polycube_name + "\n";
		log_string += "Input chartlabel: " + input_chart_label_name + "\n";
		log_string += "Input feature: " + input_feature_name + "\n";
		log_string += "file name: " + filename + "\n";
		log_string += "name prefix: " + name_prefix + "\n";
		double cut_depth = FLAGS_cd;
		int distortion_type = FLAGS_dt;
		double sigma_s = FLAGS_ss;
		double sigma_r = FLAGS_sr;
		double max_volume_ratio = FLAGS_mv;
		int cut_flag = FLAGS_cf;
		Polycube_Algorithm pa;
		TetStructure<double> input_ori, input_polycube;
		std::vector<CVec<double, 3>> ori_pts, polycube_pts;
		std::vector<unsigned int> tet_idx;
		CVec<double, 3> translation;
		input_ori.load_vtk(input_ori_name.c_str());
		std::cout << "mesh size: " << input_ori.tetra_vertices.size() << std::endl;
		input_polycube.load_vtk(input_polycube_name.c_str());
		if (input_ori.tetra_vertices.size() != input_polycube.tetra_vertices.size())
		{
			std::cout << "Input Size Error" << std::endl;
			return 1;
		}
		if (input_ori.tetra_vertices.size() == 0 || input_polycube.tetra_vertices.size() == 0)
		{
			std::cout << "Input File Empty" << std::endl;
			return 1;
		}
		pa.SetMesh(&input_ori);
		pa.SetMeshPolycube(&input_polycube);
		pa.SetSigmaR(sigma_r);
		if (!pa.load_chart_label(input_chart_label_name.c_str()))
		{
			std::cout << "Input ChartLabel Error" << std::endl;
			return 1;
		}
		if (FLAGS_input_fea_name != "nofeature.fea")
		{
			//load feature
			if (!pa.load_feature_edges(FLAGS_input_fea_name.c_str()))
			{
				std::cout << "Load Feature Error" << std::endl;
				return 1;
			}
		}
		int n_cut = -1; //0: no cut -1: input wrong -2: extraction wrong, -3: no candidate
		TetStructure<double> *pq_cut;
		
		if (cut_flag)
		{
			if (FLAGS_ct == 0)
				n_cut = pa.auto_cut_batch(cut_depth, max_volume_ratio, FLAGS_bs, FLAGS_cr, FLAGS_dt, FLAGS_ii, FLAGS_mr, FLAGS_cc, FLAGS_sh, FLAGS_cn, FLAGS_th, FLAGS_rd);
			else if (FLAGS_ct == 1)
				n_cut = pa.auto_cut_batch_min_max_depth(cut_depth, max_volume_ratio, FLAGS_bs, FLAGS_cr, FLAGS_dt, FLAGS_ii, false, FLAGS_mr, FLAGS_cc, FLAGS_sh, FLAGS_cn, FLAGS_th, FLAGS_rd);
			else if (FLAGS_ct == 2)
				n_cut = pa.auto_cut_batch_min_max_depth(cut_depth, max_volume_ratio, FLAGS_bs, FLAGS_cr, FLAGS_dt, FLAGS_ii, true, FLAGS_mr, FLAGS_cc, FLAGS_sh, FLAGS_cn, FLAGS_th, FLAGS_rd);
			if (n_cut == 0 || n_cut == -3)
			{
				pq_cut = pa.copy_mesh(pa.tet_mesh_polycube_);
				//pa.save_polycube_vtk((name_prefix + "_pqcut.vtk").c_str());
			}
			else if (n_cut == -1)
			{
				std::cout << "Input Error" << std::endl;
				return 1;
			}
			else if (n_cut == -2)
			{
				std::cout << "Extraction Error" << std::endl;
				return 1;
			}
			else
			{
				//normal case
				pq_cut = pa.copy_mesh(pa.tet_mesh_polycube_);
				//pa.save_polycube_vtk((name_prefix + "_pqcut.vtk").c_str());
				pa.auto_deformation(sigma_s, 10);
			}
			pa.save_tet_vtk((name_prefix + "_ori.vtk").c_str());
		}
		double initial_cube_length = -1.0;
		if (FLAGS_cl == 1)
		{
			std::ifstream tmp_file(input_cubelength_prefix + filename.substr(0, filename.length() - 4) + "_deform.txt"); //filename should have prefix *_cut
			if (tmp_file.is_open())
			{
				std::cout << "Cube Length File Open Successfully" << std::endl;
				tmp_file >> initial_cube_length;
				tmp_file.close();
			}
			else
			{
				//not opened 
				std::cout << "no length file found, please put it in : " << input_cubelength_prefix + filename.substr(0, filename.length() - 4) + "_deform.txt" << std::endl;
				return 1;
			}
			if (!cut_flag)
			{
				//without cut, cube size should be smaller
				initial_cube_length = initial_cube_length * 0.9;
			}
			else
			{
				initial_cube_length = initial_cube_length * 1.1;
			}
		}
		std::cout << "initial cube length: " << initial_cube_length << std::endl;
		log_string += "initial cube length: \n";
		log_string += std::to_string(initial_cube_length) + "\n";
		//bool int_solver_flag = true;
		pa.pq_flattening->set_hex_meshing_flag(FLAGS_hex_flag);
		pa.do_flattening_constrained_vtk(distortion_type, initial_cube_length, FLAGS_ld, FLAGS_ac, FLAGS_ipt, FLAGS_smd, FLAGS_is, FLAGS_fos, FLAGS_dcf);
		std::pair<int, int> flip_count = pa.compute_flips();
		//log
		log_string += "Origin Mesh Flip:\n" + std::to_string(flip_count.first) + "\n";
		log_string += "Polycube Mesh Flip:\n" + std::to_string(flip_count.second) + "\n";
		log_string += "Number of cuts:\n" + std::to_string(n_cut) + "\n";
		log_string += "Mesh Size:\n" + std::to_string(pa.tet_mesh_->tetra_vertices.size()) + " " + std::to_string(pa.tet_mesh_->tetras.size()) + "\n";
		pa.save_polycube_vtk((name_prefix + "_flattening.vtk").c_str());
		log_string += "Final Cube Length:\n" + std::to_string(pa.final_cube_length) + "\n";
		pa.save_hexex((name_prefix + "_ori.hexex").c_str(), pa.tet_mesh_, pa.tet_mesh_polycube_);
		pa.save_hexex((name_prefix + "_polycube.hexex").c_str(), pq_cut, pa.tet_mesh_polycube_);
		if (FLAGS_sf)
		{
			pa.save_flip_info((name_prefix + "_flipinfo.txt").c_str());
		}
		long end_t = clock();
		long diff1 = end_t - start_t;
		double t1 = (double)(diff1) / CLOCKS_PER_SEC;
		std::cout << "Running time: " << t1 << std::endl;
		//log
		log_string += "Running Time:\n" + std::to_string(t1) + "\n";
		log_string += "Parameters: \n";
		log_string += " ac: " + std::to_string(FLAGS_ac);
		log_string += " bs: " + std::to_string(FLAGS_bs);
		log_string += " cc: " + std::to_string(FLAGS_cc);
		log_string += " cd: " + std::to_string(FLAGS_cd);
		log_string += " cf: " + std::to_string(FLAGS_cf);
		log_string += " cl: " + std::to_string(FLAGS_cl);
		log_string += " cn: " + std::to_string(FLAGS_cn);
		log_string += " cr: " + std::to_string(FLAGS_cr);
		log_string += " ct: " + std::to_string(FLAGS_ct);
		log_string += " dff: " + std::to_string(FLAGS_dff);
		log_string += " dt: " + std::to_string(FLAGS_dt);
		log_string += " fsr: " + std::to_string(FLAGS_fsr);
		log_string += " hex_flag: " + std::to_string(FLAGS_hex_flag);
		log_string += " ii: " + std::to_string(FLAGS_ii);
		log_string += " ipt: " + std::to_string(FLAGS_ipt);
		log_string += " is: " + std::to_string(FLAGS_is);
		log_string += " ld: " + std::to_string(FLAGS_ld);
		log_string += " mr: " + std::to_string(FLAGS_mr);
		log_string += " mv: " + std::to_string(FLAGS_mv);
		log_string += " rd: " + std::to_string(FLAGS_rd);
		log_string += " sf: " + std::to_string(FLAGS_sf);
		log_string += " sh: " + std::to_string(FLAGS_sh);
		log_string += " si: " + std::to_string(FLAGS_si);
		log_string += " smd: " + std::to_string(FLAGS_smd);
		log_string += " smooth_factor: " + std::to_string(FLAGS_smooth_factor);
		log_string += " sr: " + std::to_string(FLAGS_sr);
		log_string += " ss: " + std::to_string(FLAGS_ss);
		log_string += " th: " + std::to_string(FLAGS_th);
		log_string += " offset: " + std::to_string(FLAGS_fos);
		log_string += "\n";
		std::cout << log_string;
	}
	return 1;
}
