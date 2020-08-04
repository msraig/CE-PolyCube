#include <iostream>
#include <fstream>
#include <queue>
#include <set>
#include <functional>
#include <algorithm>
#include <cassert>
#include "HexOperation.h"
#include "Hex.h"


using ig::CVec;
using namespace std;

namespace HexStructure
{
	bool read_hex_mesh_vtk(Mesh &m, char * path)
	{
		char file[300];
		sprintf(file, "%s%s", path, "");
		std::ifstream ff(file, std::ios::in);
		if (!ff.is_open()) return false;
		char s[1024], sread[1024], sread2[1024];
		Int vnum, hnum;	double x, y, z;

		bool find = false; Int lines = 0;
		while (!find)
		{
			ff.getline(s, 1023);
			if (sscanf(s, "%s %d %s", &sread, &vnum, &sread2) == 3 && (strcmp(sread, "POINTS") == 0))
				find = true;
			if (++lines > 10)
			{
				throw std::runtime_error("cannot find head of VTK!");
				return false;
			}	
		}
		/*m.V.resize(3, vnum);
		m.V.setZero();*/
		m.Vs.resize(vnum);

		//initialize vertex and cell and their relationship
		for (Int i = 0; i < vnum; i++)
		{
			ff.getline(s, 1023);
			sscanf(s, "%lf %lf %lf", &x, &y, &z);

			/*m.V(0, i) = x;
			m.V(1, i) = y;
			m.V(2, i) = z;*/

			Vertex v;
			v.v[0] = x;
			v.v[1] = y;
			v.v[2] = z;
			
			v.id = i; v.boundary = false;
			m.Vs[i] = v;
		}

		find = false;
		while (!find)
		{
			Int temp_int;
			ff.getline(s, 1023);
			if (sscanf(s, "%s %d %d", &sread, &hnum, &temp_int) == 3 && (strcmp(sread, "CELLS") == 0))
				find = true;
		}
		m.Cs.resize(hnum);
		Cell h;
		Int a, b, c, d, e, f, g, p;
		for (Int i = 0; i < hnum; i++)
			//while (ff.getline(s, 1023))
		{
			ff.getline(s, 1023);
			if (sscanf(s, "%d %d %d %d %d", &vnum, &a, &b, &c, &d) == 5 && vnum == 4) {
				h.vs.resize(4);
				//a--;b--;c--;d--;//temporarily
				h.vs[0] = a;
				h.vs[1] = b;
				h.vs[2] = c;
				h.vs[3] = d;
			}
			else if (sscanf(s, "%d %d %d %d %d %d", &vnum, &a, &b, &c, &d, &e) == 6 && vnum == 5) {
				h.vs.resize(5);
				h.vs[0] = a;
				h.vs[1] = b;
				h.vs[2] = c;
				h.vs[3] = d;
				h.vs[4] = e;
			}
			else if (sscanf(s, "%d %d %d %d %d %d %d", &vnum, &a, &b, &c, &d, &e, &f) == 7 && vnum == 6)
			{
				h.vs.resize(6);
				h.vs[0] = a;
				h.vs[1] = b;
				h.vs[2] = c;
				h.vs[3] = d;
				h.vs[4] = e;
				h.vs[5] = f;
			}
			else if (sscanf(s, "%d %d %d %d %d %d %d %d %d", &vnum, &a, &b, &c, &d, &e, &f, &g, &p) == 9 && vnum == 8)
			{
				h.vs.resize(8);
				h.vs[0] = a;
				h.vs[1] = b;
				h.vs[2] = c;
				h.vs[3] = d;
				h.vs[4] = e;
				h.vs[5] = f;
				h.vs[6] = g;
				h.vs[7] = p;
				//std::reverse(h.vs.begin(), h.vs.begin() + 4);
				//std::reverse(h.vs.begin() + 4, h.vs.end());
			}
			else {
				std::cout << "Wrong format of input file!" << endl; system("PAUSE");
			}

			h.id = i; m.Cs[h.id] = h;

			//initilize n_cs for vertex
			for (Int i = 0; i < h.vs.size(); i++) m.Vs[h.vs[i]].n_cs.push_back(h.id);
		}
		ff.close();

		return true;

	}
	
	void initial_hex_mesh(const vector<ig::CVec<double, 3>> &pts, const vector<std::vector<Int>> &idx, Mesh &m)
	{
		if (pts.size() == 0 || idx.size() == 0)
		{
			cout << "Empty input array" << endl;
			return;
		}

		Int vnum((Int)pts.size()), hnum((Int)idx.size());

		m.Vs.resize(vnum);

		for (Int i = 0; i < vnum; i++)
		{

			Vertex v;
			v.v[0] = pts[i][0];
			v.v[1] = pts[i][1];
			v.v[2] = pts[i][2];

			v.id = i; v.boundary = false;
			m.Vs[i] = v;
		}

		m.Cs.resize(hnum);
		Cell h;
		for (Int i = 0; i < hnum; i++)
			//while (ff.getline(s, 1023))
		{
			{
				h.vs = idx[i];
			}

			h.id = i; m.Cs[h.id] = h;

			//initilize n_cs for vertex
			for (Int i = 0; i < h.vs.size(); i++) m.Vs[h.vs[i]].n_cs.push_back(h.id);
		}



	}

	void get_hex_vert_cell(const Mesh &m, vector<ig::CVec<Real, 3>> &pts, vector<std::vector<Int>> &idx)
	{
		pts.clear();
		idx.clear();
		for (Int i = 0; i < m.Vs.size(); i++)
		{
			pts.push_back(m.Vs[i].v);
		}
		for (Int i = 0; i < m.Cs.size(); i++)
		{
			idx.push_back(m.Cs[i].vs);
		}
	}

	void write_hex_mesh_vtk(const Mesh &m, char * path)
	{
		std::fstream f(path, std::ios::out);

		f << "# vtk DataFile Version 2.0" << std::endl << "mesh vtk data - converted from .off" << std::endl;
		f << "ASCII" << std::endl;
		f << "DATASET UNSTRUCTURED_GRID" << std::endl;

		f << "POINTS " << m.Vs.size() << " double" << std::endl;
		for (Int i = 0; i < m.Vs.size(); i++)
			//f << m.V(0, i) << " " << m.V(1, i) << " " << m.V(2, i) << std::endl;
			f << m.Vs[i].v[0] << " " << m.Vs[i].v[1] << " " << m.Vs[i].v[2] << std::endl;

		{
			size_t vnum = m.Cs[0].vs.size();
			f << "CELLS " << m.Cs.size() << " " << m.Cs.size() * (vnum + 1) << std::endl;

			for (Int i = 0; i < m.Cs.size(); i++)
			{
				f << vnum << " ";
				for (Int j = 0; j < vnum; j++)
					f << m.Cs[i].vs[j] << " ";
				/*f << hmi.Hs[i].vs[3] << " ";
				f << hmi.Hs[i].vs[2] << " ";
				f << hmi.Hs[i].vs[1] << " ";
				f << hmi.Hs[i].vs[0] << " ";
				f << hmi.Hs[i].vs[7] << " ";
				f << hmi.Hs[i].vs[6] << " ";
				f << hmi.Hs[i].vs[5] << " ";
				f << hmi.Hs[i].vs[4] << " ";*/
				f << std::endl;
			}
			f << "CELL_TYPES " << m.Cs.size() << std::endl;
			for (Int i = 0; i < m.Cs.size(); i++)
					f << 12 << std::endl;
		}


		/*f << "POINT_DATA " << m.Vs.size() << std::endl;
		f << "SCALARS fixed int" << std::endl;
		f << "LOOKUP_TABLE default" << std::endl;
		for (Int i = 0; i < m.Vs.size(); i++) {
			if (m.Vs[i].boundary) f << "1" << std::endl; else f << "0" << std::endl;
		}*/
		f.close();
	}

	void write_singularity_vtk(const Mesh &mesh, const Singularity &si, char * path)
	{
		std::fstream ff(path, std::ios::out);
		ff << "# vtk DataFile Version 2.0" << std::endl << "mesh vtk data" << std::endl;
		ff << "ASCII" << std::endl;
		ff << "DATASET POLYDATA" << std::endl;


		ff << "POINTS " << mesh.Vs.size() << " double" << std::endl;
		for (Int i = 0; i < mesh.Vs.size(); i++)
			//ff << mesh.V(0, i) << " " << mesh.V(1, i) << " " << mesh.V(2, i) << std::endl;
			ff << mesh.Vs[i].v[0] << " " << mesh.Vs[i].v[1] << " " << mesh.Vs[i].v[2] << std::endl;

		ff << "VERTICES " << si.SVs.size() << " " << si.SVs.size() * 2 << endl;
		for (Int i = 0; i < si.SVs.size(); i++) ff << "1 " << si.SVs[i].oid << endl;

		size_t line_num = 0;
		for (Int i = 0; i < si.SEs.size(); i++) line_num += si.SEs[i].vs_onmesh.size();

		ff << "LINES " << si.SEs.size() << " " << line_num + si.SEs.size() << endl;
		for (Int i = 0; i < si.SEs.size(); i++) {
			ff << si.SEs[i].vs_onmesh.size() << " ";
			for (Int j = 0; j < si.SEs[i].vs_onmesh.size(); j++) ff << si.SEs[i].vs_onmesh[j] << " ";
			ff << endl;
		}


		ff << "POINT_DATA " << mesh.Vs.size() << std::endl;
		ff << "SCALARS V_Scalars int" << std::endl;
		ff << "LOOKUP_TABLE V_Table" << std::endl;
		std::vector<short> V_tag(mesh.Vs.size(), 0);
		for (Int i = 0; i < mesh.Vs.size(); i++) if (mesh.Vs[i].boundary) V_tag[i] = 1;
		for (Int j = 0; j < si.SVs.size(); j++) V_tag[si.SVs[j].oid] = 3;
		for (Int j = 0; j < si.SEs.size(); j++) for (Int k = 0; k < si.SEs[j].vs_onmesh.size(); k++)
			V_tag[si.SEs[j].vs_onmesh[k]] = 3;

		for (Int i = 0; i < V_tag.size(); i++) ff << V_tag[i] << std::endl;

		ff.close();
	}
	void write_metamesh_vtk(const Mesh &mesh, const MetaMesh &mm, char * path)
	{
		//get meta-mesh from metamesh, face may not be coincident
		std::fstream ff(path, std::ios::out);
		ff << "# vtk DataFile Version 2.0" << std::endl << "mesh vtk data - converted from .off" << std::endl;
		ff << "ASCII" << std::endl;
		ff << "DATASET UNSTRUCTURED_GRID" << std::endl;


		ff << "POINTS " << mm.FVs.size() << " double" << std::endl;
		for (Int i = 0; i < mm.FVs.size(); i++)
			//ff << mesh.V(0, mm.FVs[i].hid) << " " << mesh.V(1, mm.FVs[i].hid) << " " << mesh.V(2, mm.FVs[i].hid) << std::endl;
			ff << mesh.Vs[mm.FVs[i].hid].v[0] << " " << mesh.Vs[mm.FVs[i].hid].v[1] << " " << mesh.Vs[mm.FVs[i].hid].v[2] << std::endl;


		//output hex below
		size_t vnum = mm.FHs[0].vs.size();
		ff << "CELLS " << mm.FHs.size() << " " << mm.FHs.size() * (vnum + 1) << std::endl;

		size_t n_hex = mm.FHs.size();
		for (size_t i = 0; i < n_hex; i++)
		{
			ff << vnum << " ";
			for (size_t j = 0; j < vnum; j++)
			{
				ff << mm.FHs[i].vs[j] << " ";
			}
			ff << std::endl;
		}

		ff << "CELL_TYPES " << n_hex << std::endl;
		for (size_t i = 0; i < n_hex; i++)
		{
			ff << 12 << std::endl;
		}

		ff.close();
	}
	

	void build_connectivity(Mesh &m)
	{
		m.Es.clear(); if (m.Cs.size()) m.Fs.clear();

		std::vector<std::vector<Int>> total_fs(m.Cs.size() * 6);
		std::vector<std::tuple<Int, Int, Int, Int, Int, Int, Int>> tempF(m.Cs.size() * 6);
		std::vector<Int> vs(4);
		for (Int i = 0; i < m.Cs.size(); ++i) {
			for (short j = 0; j < 6; j++) {
				for (short k = 0; k < 4; k++) vs[k] = m.Cs[i].vs[hex_face_table[j][k]];
				Int id = 6 * i + j;
				total_fs[id] = vs;
				std::sort(vs.begin(), vs.end());
				tempF[id] = std::make_tuple(vs[0], vs[1], vs[2], vs[3], id, i, j);
			}
			m.Cs[i].fs.resize(6);
		}
		//key part here, comparison is meaningful only when sorting is done
		std::sort(tempF.begin(), tempF.end());
		m.Fs.reserve(tempF.size() / 3);
		Face f; f.boundary = true;
		Int F_num = 0;
		for (Int i = 0; i < tempF.size(); ++i) {
			if (i == 0 || (i != 0 &&
				(std::get<0>(tempF[i]) != std::get<0>(tempF[i - 1]) || std::get<1>(tempF[i]) != std::get<1>(tempF[i - 1]) ||
					std::get<2>(tempF[i]) != std::get<2>(tempF[i - 1]) || std::get<3>(tempF[i]) != std::get<3>(tempF[i - 1])))) {
				//first face || (at least one vertex df from the former)
				f.id = F_num; F_num++;
				f.vs = total_fs[std::get<4>(tempF[i])];
				m.Fs.push_back(f);
			}
			else if (i != 0 && (std::get<0>(tempF[i]) == std::get<0>(tempF[i - 1]) && std::get<1>(tempF[i]) == std::get<1>(tempF[i - 1]) &&
				std::get<2>(tempF[i]) == std::get<2>(tempF[i - 1]) && std::get<3>(tempF[i]) == std::get<3>(tempF[i - 1])))
				m.Fs[F_num - 1].boundary = false;  //duplicate face found, no boundary

			m.Cs[std::get<5>(tempF[i])].fs[std::get<6>(tempF[i])] = F_num - 1;
		}

		std::vector<std::tuple<Int, Int, Int, Int>> temp(m.Fs.size() * 4);
		for (Int i = 0; i < m.Fs.size(); ++i) {
			for (Int j = 0; j < 4; ++j) {
				Int v0 = m.Fs[i].vs[j], v1 = m.Fs[i].vs[(j + 1) % 4];
				if (v0 > v1) std::swap(v0, v1);
				temp[4 * i + j] = std::make_tuple(v0, v1, i, j);
			}
			m.Fs[i].es.resize(4);
		}
		std::sort(temp.begin(), temp.end());
		m.Es.reserve(temp.size() / 2);
		Int E_num = 0;
		Edge e; e.boundary = false; e.vs.resize(2);
		for (Int i = 0; i < temp.size(); ++i) {
			if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) ||
				std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
				e.id = E_num; E_num++;
				e.vs[0] = std::get<0>(temp[i]);
				e.vs[1] = std::get<1>(temp[i]);
				m.Es.push_back(e);
			}
			m.Fs[std::get<2>(temp[i])].es[std::get<3>(temp[i])] = E_num - 1;
		}
		//boundary
		for (auto &v : m.Vs) v.boundary = false;
		for (Int i = 0; i < m.Fs.size(); ++i)
			if (m.Fs[i].boundary) for (Int j = 0; j < 4; ++j) {
				Int eid = m.Fs[i].es[j];
				m.Es[eid].boundary = true;
				m.Vs[m.Es[eid].vs[0]].boundary = m.Vs[m.Es[eid].vs[1]].boundary = true;
			}
		
		//neighboring cell of face
		for (Int i = 0; i < m.Cs.size(); i++) {
			for (Int j = 0; j < m.Cs[i].fs.size(); j++) m.Fs[m.Cs[i].fs[j]].n_cs.push_back(i);
		}

		//neighboring face of edge and vertex
		for (Int i = 0; i < m.Fs.size(); i++) {
			for (Int j = 0; j < m.Fs[i].es.size(); j++) m.Es[m.Fs[i].es[j]].n_fs.push_back(i);
			for (Int j = 0; j < m.Fs[i].vs.size(); j++) m.Vs[m.Fs[i].vs[j]].n_fs.push_back(i);
		}

		//neighboring edge and neighboring vertex of vertex
		for (Int i = 0; i < m.Es.size(); i++) {
			Int v0 = m.Es[i].vs[0], v1 = m.Es[i].vs[1];
			m.Vs[v0].n_es.push_back(i);
			m.Vs[v1].n_es.push_back(i);
			m.Vs[v0].n_vs.push_back(v1);
			m.Vs[v1].n_vs.push_back(v0);
		}
		
		//neighboring cells of edge
		for (Int i = 0; i < m.Es.size(); i++) {
			std::vector<Int> nhs;
			for (Int j = 0; j < m.Es[i].n_fs.size(); j++) {
				Int nfid = m.Es[i].n_fs[j];
				nhs.insert(nhs.end(), m.Fs[nfid].n_cs.begin(), m.Fs[nfid].n_cs.end());
			}
			std::sort(nhs.begin(), nhs.end()); nhs.erase(std::unique(nhs.begin(), nhs.end()), nhs.end());
			m.Es[i].n_cs = nhs;
		}
	}
	void extract_singularity(Mesh &mesh, Singularity &si)
	{
		si.SVs.clear(); si.SEs.clear();

		Int INVALID_V = (Int)-1;
		Int INVALID_E = (Int)-1;
		std::vector<Int> V_flag(mesh.Vs.size(), INVALID_V), E_flag(mesh.Es.size(), false);

		for (auto &v : mesh.Vs) { v.mvid = -1; v.svid = -1; }

		Int SV_count = 0, SE_count = 0;
		vector<Singularity_E> circle_ses;//circular singular edges

		for (int i = 0; i < mesh.Es.size(); i++) {
			if (E_flag[i]) continue;

			if ((!mesh.Es[i].boundary && mesh.Es[i].n_cs.size() == Interior_RegularE) ||
				(mesh.Es[i].boundary && mesh.Es[i].n_cs.size() == Boundary_RegularE))
				continue;//non-singular e-->continue;

			std::function<bool(Int, Int, Int &)> singular_proceed = [&](Int vid, Int eid, Int &neid)->bool {
				Int num1 = 0, num2 = 0;
				for (Int j = 0; j < mesh.Vs[vid].n_es.size(); j++) {
					Int cur_e = mesh.Vs[vid].n_es[j];
					if (cur_e == eid) continue;
					if ((!mesh.Es[cur_e].boundary && mesh.Es[cur_e].n_cs.size() != Interior_RegularE) ||
						(mesh.Es[cur_e].boundary && mesh.Es[cur_e].n_cs.size() != Boundary_RegularE))
						num1++;

					if (mesh.Es[cur_e].boundary == mesh.Es[eid].boundary&&mesh.Es[cur_e].n_cs.size() == mesh.Es[eid].n_cs.size()) {
						num2++; neid = cur_e;
					}
				}
				if (num1 == 1 && num2 == 1)
					return true;
				return false;
			};
			Int v_left = mesh.Es[i].vs[0], v_right = mesh.Es[i].vs[1];
			Int sv_left, sv_right;
			std::vector<Int> vs_left, vs_right, es_left, es_right;

			bool is_circle = false;
			//left 
			es_left.push_back(i); vs_left.push_back(v_left);
			Int cur_e = i, next_e = INVALID_E;
			while (singular_proceed(v_left, cur_e, next_e)) {
				cur_e = next_e;
				if (cur_e == i) { is_circle = true; break; }
				es_left.push_back(next_e);
				if (mesh.Es[cur_e].vs[0] == v_left) v_left = mesh.Es[cur_e].vs[1]; else v_left = mesh.Es[cur_e].vs[0];
				vs_left.push_back(v_left);
			}

			if (is_circle) {
				Singularity_E se; se.id = (Int)si.SEs.size();
				se.circle = true;
				se.es_onmesh = es_left; se.vs_onmesh = vs_left;
				for (Int j = 0; j < se.es_onmesh.size(); j++) E_flag[se.es_onmesh[j]] = true;
				se.boundary = mesh.Es[i].boundary;
				si.SEs.push_back(se);
				circle_ses.push_back(se);
				continue;
			}

			sv_left = V_flag[v_left];
			if (V_flag[v_left] == INVALID_V) {
				sv_left = (Int)si.SVs.size(); V_flag[v_left] = sv_left;
				Singularity_V sv; sv.fake = false;
				sv.id = sv_left;
				sv.oid = v_left;
				sv.boundary = mesh.Vs[v_left].boundary;
				si.SVs.push_back(sv);
			}
			//right
			vs_right.push_back(v_right);
			cur_e = i, next_e = INVALID_E;
			while (singular_proceed(v_right, cur_e, next_e)) {
				cur_e = next_e;
				if (mesh.Es[cur_e].vs[0] == v_right) v_right = mesh.Es[cur_e].vs[1]; else v_right = mesh.Es[cur_e].vs[0];
				vs_right.push_back(v_right);
				es_right.push_back(next_e);
			}
			sv_right = V_flag[v_right];
			if (V_flag[v_right] == INVALID_V) {
				sv_right = (Int)si.SVs.size(); V_flag[v_right] = sv_right;
				Singularity_V sv; sv.fake = false;
				sv.id = sv_right;
				sv.oid = v_right;
				sv.boundary = mesh.Vs[v_right].boundary;
				si.SVs.push_back(sv);
			}
			//se
			Singularity_E se; se.id = (Int)si.SEs.size(); se.circle = false;
			std::reverse(vs_left.begin(), vs_left.end());
			se.vs_onmesh = vs_left; se.vs_onmesh.insert(se.vs_onmesh.end(), vs_right.begin(), vs_right.end());
			std::reverse(es_left.begin(), es_left.end());
			se.es_onmesh = es_left; se.es_onmesh.insert(se.es_onmesh.end(), es_right.begin(), es_right.end());
			for (Int j = 0; j < se.es_onmesh.size(); j++) E_flag[se.es_onmesh[j]] = true;
			se.boundary = mesh.Es[i].boundary;

			if (sv_left == sv_right) {
				se.vs.resize(1); se.vs[0] = sv_left;
				se.vs_onmesh.erase(se.vs_onmesh.begin() + se.vs_onmesh.size() - 1);
				se.circle = true;
				circle_ses.push_back(se);
			}
			else {
				se.vs.resize(2); se.vs[0] = sv_left; se.vs[1] = sv_right;
			}
			si.SEs.push_back(se);
		}

		for (auto sv : si.SVs) mesh.Vs[sv.oid].svid = sv.id;
	}

	//Preparation for extraction of metamesh
	enum V_tag {
		R = -5,//regular
		S,//singular node
		E,//extra-ordinary: on singular edge
		Bn,//base-complex node
		Be//on base-complex edge
	};

	void node_on_singularity_circle(Singularity &si, MetaMesh &mm, Mesh &mesh, vector<Int> &Nodes) {

		//nodes on circular singularities
		std::vector<int> v_flag(mesh.Vs.size(), V_tag::R);
		std::vector<Int> circulars;
		for (Int i = 0; i < si.SEs.size(); i++) {
			if (si.SEs[i].circle) circulars.push_back(i);
			for (Int j = 0; j < si.SEs[i].vs_onmesh.size(); j++)
				if (si.SEs[i].circle) v_flag[si.SEs[i].vs_onmesh[j]] = V_tag::E;
				else v_flag[si.SEs[i].vs_onmesh[j]] = i;
		}
		for (Int i = 0; i < si.SVs.size(); i++) v_flag[si.SVs[i].oid] = V_tag::S;

		vector<bool> E_tag(mesh.Es.size(), false);

		vector<vector<Int>> Node_cans(circulars.size());

		function<bool()>intersect_singularVs = [&]() -> bool {
			bool new_S = false;

			for (Int i = 0; i < circulars.size(); i++) {

				if (Node_cans[i].size() >= 3)continue;

				auto &vs_link = si.SEs[circulars[i]].vs_onmesh;
				auto &es_link = si.SEs[circulars[i]].es_onmesh;

				for (auto eid : es_link)E_tag[eid] = true;

				for (Int j = 0; j < vs_link.size(); j++) {
					auto vid = vs_link[j];
					if (v_flag[vid] == V_tag::S) continue;

					for (auto eid : mesh.Vs[vid].n_es) {
						if (E_tag[eid]) continue;

						Int vid_pre = vid, vid_aft = mesh.Es[eid].vs[0]; if (vid_aft == vid_pre) vid_aft = mesh.Es[eid].vs[1];
						bool found = false;

						while (true) {
							if (v_flag[vid_aft] == V_tag::S) {
								Node_cans[i].push_back(vid); v_flag[vid] = V_tag::S;
								found = true; new_S = true; break;
							}

							bool find_next = false;
							for (auto nvid : mesh.Vs[vid_aft].n_vs) {
								if (nvid == vid_pre || nvid == vid || v_flag[nvid] == V_tag::E || v_flag[nvid] >= 0) continue;
								std::sort(mesh.Vs[nvid].n_fs.begin(), mesh.Vs[nvid].n_fs.end());
								std::sort(mesh.Vs[vid_pre].n_fs.begin(), mesh.Vs[vid_pre].n_fs.end());
								vector<Int> common_fs;
								std::set_intersection(mesh.Vs[vid_pre].n_fs.begin(), mesh.Vs[vid_pre].n_fs.end(),
									mesh.Vs[nvid].n_fs.begin(), mesh.Vs[nvid].n_fs.end(), std::back_inserter(common_fs));
								if (!common_fs.size()) {
									vid_pre = vid_aft; vid_aft = nvid; find_next = true; break;
								}
							}
							if (!find_next) break;
						}
						if (found) break;
					}
				}
				for (auto eid : es_link)E_tag[eid] = false;
			}
			return new_S;
		};
		function<bool()>intersect_singularEVs = [&]() -> bool {
			bool new_S = false;

			for (Int i = 0; i < circulars.size(); i++) {

				if (Node_cans[i].size() >= 3)continue;

				auto &vs_link = si.SEs[circulars[i]].vs_onmesh;
				auto &es_link = si.SEs[circulars[i]].es_onmesh;

				for (auto eid : es_link)E_tag[eid] = true;

				vector<pair<Int, Int>> candidates; vector<int32_t> se_vote(si.SEs.size(), 0);

				for (Int j = 0; j < vs_link.size(); j++) {
					auto vid = vs_link[j];
					if (v_flag[vid] == V_tag::S) continue;

					for (auto eid : mesh.Vs[vid].n_es) {
						if (E_tag[eid]) continue;

						Int vid_pre = vid, vid_aft = mesh.Es[eid].vs[0]; if (vid_aft == vid_pre) vid_aft = mesh.Es[eid].vs[1];
						bool found = false;

						while (true) {
							if (v_flag[vid_aft] >= 0) {
								candidates.push_back(make_pair(vid, v_flag[vid_aft]));
								se_vote[v_flag[vid_aft]]++;
								break;
							}

							bool find_next = false;
							for (auto nvid : mesh.Vs[vid_aft].n_vs) {
								if (nvid == vid_pre || nvid == vid || v_flag[nvid] == V_tag::E) continue;
								std::sort(mesh.Vs[nvid].n_fs.begin(), mesh.Vs[nvid].n_fs.end());
								std::sort(mesh.Vs[vid_pre].n_fs.begin(), mesh.Vs[vid_pre].n_fs.end());
								vector<Int> common_fs;
								std::set_intersection(mesh.Vs[vid_pre].n_fs.begin(), mesh.Vs[vid_pre].n_fs.end(),
									mesh.Vs[nvid].n_fs.begin(), mesh.Vs[nvid].n_fs.end(), std::back_inserter(common_fs));
								if (!common_fs.size()) {
									vid_pre = vid_aft; vid_aft = nvid; find_next = true; break;
								}
							}
							if (!find_next) break;
						}
						if (found) break;
					}
				}
				for (auto eid : es_link)E_tag[eid] = false;

				for (Int j = 0; j < se_vote.size(); j++) if (se_vote[j] > 0 && se_vote[j] <= 2) {//if (se_vote[j] > 0 && se_vote[j] <= si.SEs[j].vs_link.size()) {//
					for (auto a_pair : candidates) if (a_pair.second == j) {
						Node_cans[i].push_back(a_pair.first); v_flag[a_pair.first] = V_tag::S;
						new_S = true;
					}
				}
				if (new_S) {
					std::sort(Node_cans[i].begin(), Node_cans[i].end());
					Node_cans[i].erase(std::unique(Node_cans[i].begin(), Node_cans[i].end()), Node_cans[i].end());
					return true;
				}
			}
			return false;
		};
		function<bool()>intersect_parallelVs = [&]() -> bool {
			bool new_S = false;

			for (Int i = 0; i < circulars.size(); i++) {

				if (Node_cans[i].size() >= 3)continue;

				auto &vs_link = si.SEs[circulars[i]].vs_onmesh;
				auto &es_link = si.SEs[circulars[i]].es_onmesh;

				size_t len = vs_link.size() / 3;
				vector<Int> candidates(3);
				candidates[0] = vs_link[0];
				candidates[1] = vs_link[len];
				candidates[2] = vs_link[2 * len];

				for (Int j = 0; j < 3; j++) {
					if (v_flag[candidates[j]] == V_tag::S) continue;
					Node_cans[i].push_back(candidates[j]);
					v_flag[candidates[j]] = V_tag::S;
					new_S = true;
					if (Node_cans[i].size() >= 3) break;
				}
				return new_S;
			}
			return false;
		};

		for (Int i = 0; i < circulars.size(); i++) {

			auto &vs_link = si.SEs[circulars[i]].vs_onmesh;
			auto &es_link = si.SEs[circulars[i]].es_onmesh;

			for (Int j = 0; j < vs_link.size(); j++) {
				auto vid = vs_link[j];
				if (v_flag[vid] == V_tag::S) Node_cans[i].push_back(vid);
			}
		}

		bool still_circular = true;
		while (still_circular) {
			while (intersect_singularVs());
			while (intersect_singularEVs()) {
				while (intersect_singularVs());
			}
			while (intersect_parallelVs()) {
				while (intersect_singularVs());
				while (intersect_singularVs());
				while (intersect_singularEVs()) {
					while (intersect_singularVs());
				}
			}
			still_circular = false;
			for (auto can : Node_cans)if (can.size() <= 2) {
				still_circular = true; break;
			}
		}

		std::fill(v_flag.begin(), v_flag.end(), V_tag::R);
		for (Int i = 0; i < si.SVs.size(); i++) v_flag[si.SVs[i].oid] = V_tag::S;
		Nodes.clear();
		for (auto cans : Node_cans) for (auto vid : cans) if (v_flag[vid] != V_tag::S)  Nodes.push_back(vid);
		std::sort(Nodes.begin(), Nodes.end());
		Nodes.erase(std::unique(Nodes.begin(), Nodes.end()), Nodes.end());
	}

	void node_edge_extraction(Singularity &si, MetaMesh &mm, Mesh &mesh) {
		size_t NON_SE = mesh.Vs.size() + 1, MULTIPLE_SE = mesh.Vs.size(), INVALIDE_FE = mesh.Es.size();

		mm.FVs.clear(); mm.FEs.clear();

		std::vector<V_tag> v_tag(mesh.Vs.size(), V_tag::R);
		std::vector<size_t> v_neibor_se(mesh.Vs.size(), NON_SE), v_on_fe(mesh.Vs.size(), INVALIDE_FE);
		std::vector<bool> e_flag(mesh.Es.size(), false);

		for (Int i = 0; i < si.SEs.size(); i++)
			for (Int j = 0; j < si.SEs[i].vs_onmesh.size(); j++) {
				Int vid = si.SEs[i].vs_onmesh[j]; v_tag[vid] = V_tag::E;
				if (v_neibor_se[vid] != NON_SE) v_neibor_se[vid] = MULTIPLE_SE;
				else v_neibor_se[vid] = i;
			}
		//nodes
		vector<Int> nodes(si.SVs.size());
		for (Int i = 0; i < si.SVs.size(); i++) nodes[i] = si.SVs[i].oid;
		//nodes on circular singularities
		vector<Int> Snodes;
		node_on_singularity_circle(si, mm, mesh, Snodes);
		if (Snodes.size()) nodes.insert(nodes.end(), Snodes.begin(), Snodes.end());
		//node_pool
		std::queue<Int> node_pool;
		for (Int i = 0; i < nodes.size(); i++) {
			MetaMesh_V fv;
			fv.id = (Int)mm.FVs.size(); fv.hid = nodes[i];
			v_tag[fv.hid] = V_tag::Bn;
			mesh.Vs[fv.hid].mvid = fv.id; mm.FVs.push_back(fv);
			node_pool.push(fv.id);
		}
		//return;
		//nodes --> edges
		while (!node_pool.empty()) {
			Int fvid = node_pool.front(); node_pool.pop();
			Int vid = mm.FVs[fvid].hid;
			for (Int i = 0; i < mesh.Vs[vid].n_es.size(); i++) {
				Int eid = mesh.Vs[vid].n_es[i];
				if (e_flag[eid]) continue; e_flag[eid] = true;

				Int vid_next = mesh.Es[eid].vs[0]; if (vid_next == vid) vid_next = mesh.Es[eid].vs[1];

				MetaMesh_E fe; fe.id = (Int)mm.FEs.size();
				fe.vs.push_back(fvid);
				fe.vs_onmesh.push_back(vid); fe.vs_onmesh.push_back(vid_next);
				fe.es_onmesh.push_back(eid);
				//from fvid, trace along eid direction...
				Int vid_pre = vid;
				while (true) {
					if (v_tag[vid_next] == V_tag::Bn) {//hit another Bc node
						fe.vs.push_back(mesh.Vs[vid_next].mvid);
						mm.FEs.push_back(fe);
						break;
					}
					else if (v_tag[vid_next] == V_tag::Be) {//hit a vertex on Be
															//new node
						MetaMesh_V fv; fv.id = (Int)mm.FVs.size(); fv.hid = vid_next;
						v_tag[fv.hid] = V_tag::Bn;
						mesh.Vs[fv.hid].mvid = fv.id; mm.FVs.push_back(fv);
						node_pool.push(fv.id);
						//edge end
						fe.vs.push_back(fv.id); mm.FEs.push_back(fe);
						//split feid --> feid, feid1
						Int feid = (Int)v_on_fe[vid_next];
						Int feid_v1 = mm.FEs[feid].vs[1]; mm.FEs[feid].vs[1] = fv.id;

						MetaMesh_E fe_new; fe_new.id = (Int)mm.FEs.size();
						fe_new.vs.push_back(fv.id); fe_new.vs.push_back(feid_v1);

						std::vector<Int> &vs_link = mm.FEs[feid].vs_onmesh, &es_link = mm.FEs[feid].es_onmesh;
						size_t vid_next_pos = (size_t)(std::find(vs_link.begin(), vs_link.end(), vid_next) - vs_link.begin());
						fe_new.vs_onmesh.insert(fe_new.vs_onmesh.end(), vs_link.begin() + vid_next_pos, vs_link.end());
						fe_new.es_onmesh.insert(fe_new.es_onmesh.end(), es_link.begin() + vid_next_pos, es_link.end());
						vs_link.erase(vs_link.begin() + vid_next_pos + 1, vs_link.end());
						es_link.erase(es_link.begin() + vid_next_pos, es_link.end());

						for (Int j = 0; j < fe_new.vs_onmesh.size(); j++) v_on_fe[fe_new.vs_onmesh[j]] = fe_new.id;
						mm.FEs.push_back(fe_new);
						break;
					}
					else if (v_tag[vid_next] == V_tag::E && (v_neibor_se[vid_pre] != MULTIPLE_SE && v_neibor_se[vid_next] != v_neibor_se[vid_pre])) {//hit a vertex on Se
																											//new node
						MetaMesh_V fv; fv.id = (Int)mm.FVs.size(); fv.hid = vid_next;
						v_tag[fv.hid] = V_tag::Bn;
						mesh.Vs[fv.hid].mvid = fv.id; mm.FVs.push_back(fv);
						node_pool.push(fv.id);
						//edge end
						fe.vs.push_back(fv.id); mm.FEs.push_back(fe);
						break;
					}

					v_tag[vid_next] = V_tag::Be;
					v_on_fe[vid_next] = fe.id;
					bool find_next = false;
					for (Int j = 0; j < mesh.Vs[vid_next].n_es.size(); j++) {
						Int neid = mesh.Vs[vid_next].n_es[j];
						Int nvid = mesh.Es[neid].vs[0]; if (nvid == vid_next) nvid = mesh.Es[neid].vs[1];
						if (nvid == vid_pre) continue;

						bool found = true;
						for (Int k = 0; k < mesh.Es[neid].n_cs.size(); k++) {
							Int nhid = mesh.Es[neid].n_cs[k];
							if (std::find(mesh.Cs[nhid].vs.begin(), mesh.Cs[nhid].vs.end(), vid_pre) != mesh.Cs[nhid].vs.end()) {
								found = false; break;
							}
						}
						if (found) {
							vid_pre = vid_next; vid_next = nvid;
							fe.vs_onmesh.push_back(vid_next);
							fe.es_onmesh.push_back(neid);
							find_next = true;
							e_flag[neid] = true;
							break;
						}
					}
					if (!find_next) {
						//new node
						MetaMesh_V fv; fv.id = (Int)mm.FVs.size(); fv.hid = vid_next;
						v_tag[fv.hid] = V_tag::Bn;
						mesh.Vs[fv.hid].mvid = fv.id; mm.FVs.push_back(fv);
						node_pool.push(fv.id);
						//edge end
						fe.vs.push_back(fv.id); mm.FEs.push_back(fe);
						break;
					}
				}
			}
		}

		for (Int i = 0; i < mm.FEs.size(); i++) {
			Int v0 = mm.FEs[i].vs[0], v1 = mm.FEs[i].vs[1];
			mm.FVs[v0].n_fvs.push_back(v1);
			mm.FVs[v1].n_fvs.push_back(v0);
			mm.FVs[v0].n_fes.push_back(i);
			mm.FVs[v1].n_fes.push_back(i);
		}

		for (auto fv : mm.FVs)mesh.Vs[fv.hid].mvid = fv.id;

		//char path[300] = "../datasets/orig.hex_frame.vtk";
		//h_io io;
		//io.write_Frame_VTK(frame,mesh,path);
	}

	bool face_extraction(Singularity &si, MetaMesh &mm, Mesh &mesh)
	{
		mm.FFs.clear();
		size_t INVALID_E = mm.FEs.size(), INVALID_F = mesh.Fs.size();
		std::vector<size_t> e_tag(mesh.Es.size(), INVALID_E);
		for (Int i = 0; i < mm.FEs.size(); i++)
			for (Int j = 0; j < mm.FEs[i].es_onmesh.size(); j++) e_tag[mm.FEs[i].es_onmesh[j]] = i;
		std::vector<size_t> f_flag(mesh.Fs.size(), INVALID_F);
		std::vector<bool> fe_flag(mm.FEs.size(), false);

		//edge info should be given
		for (Int i = 0; i < mm.FEs.size(); i++) {
			Int eid = mm.FEs[i].es_onmesh[0]; // first edge of the metamesh edge
			for (Int j = 0; j < mesh.Es[eid].n_fs.size(); j++) {
				Int fid = mesh.Es[eid].n_fs[j]; //one of the neighbor faces of edge eid
				if (f_flag[fid] != INVALID_F) continue; //only consider f which havn't been processed

				MetaMesh_F ff; ff.id = (Int)mm.FFs.size();
				ff.boundary = mesh.Fs[fid].boundary;
				ff.es.push_back(i); fe_flag[i] = true;

				std::queue<Int> f_pool; f_pool.push(fid); //push the first face in the queue
				while (!f_pool.empty()) {
					fid = f_pool.front(); f_pool.pop();
					if (f_flag[fid] != INVALID_F) continue;
					f_flag[fid] = ff.id;
					ff.ffs_grid.push_back(fid);

					for (Int k = 0; k < 4; k++) {
						Int feid = mesh.Fs[fid].es[k];
						if (e_tag[feid] != INVALID_E) {
							//edge have be processed
							if (!fe_flag[e_tag[feid]]) {
								fe_flag[e_tag[feid]] = true;
								ff.es.push_back((Int)e_tag[feid]);
							}
							continue;
						}
						for (Int m = 0; m < mesh.Es[feid].n_fs.size(); m++) {
							Int enfid = mesh.Es[feid].n_fs[m];
							if (f_flag[enfid] != INVALID_F) continue;
							bool pass = true;
							for (Int n = 0; n < mesh.Fs[enfid].n_cs.size(); n++) {
								Int fnhid = mesh.Fs[enfid].n_cs[n];
								for (Int p = 0; p < 6; p++) {
									Int cur_fid = mesh.Cs[fnhid].fs[p];
									if (f_flag[cur_fid] == ff.id) { pass = false; break; }
								}
								if (!pass) break;
							}
							if (pass) f_pool.push(enfid);
						}
					}
				}
				mm.FFs.push_back(ff);
				for (Int k = 0; k < ff.es.size(); k++) fe_flag[ff.es[k]] = false;
			}
		}
		//re-order es
		for (Int i = 0; i < mm.FFs.size(); i++) {

			if (mm.FFs[i].es.size() != 4) return false;

			Int e0 = mm.FFs[i].es[0], e1 = -1, e2 = -1, e3 = -1;
			Int v0 = mm.FEs[e0].vs[0], v1 = mm.FEs[e0].vs[1], v2 = -1, v3 = -1;
			//e1, v2
			for (Int j = 1; j < 4; j++) {
				Int eid = mm.FFs[i].es[j], v0_ = mm.FEs[eid].vs[0], v1_ = mm.FEs[eid].vs[1];
				if (v0_ == v1 || v1_ == v1) {
					if (v0_ == v1) v2 = v1_; else v2 = v0_;
					e1 = eid; break;
				}
			}
			//e3, v3
			for (Int j = 1; j < 4; j++) {
				Int eid = mm.FFs[i].es[j], v0_ = mm.FEs[eid].vs[0], v1_ = mm.FEs[eid].vs[1];
				if (v0_ == v0 || v1_ == v0) {
					if (v0_ == v0) v3 = v1_; else v3 = v0_;
					e3 = eid; break;
				}
			}


			for (Int j = 1; j < 4; j++)
				if (mm.FFs[i].es[j] != e1 && mm.FFs[i].es[j] != e3)
					e2 = mm.FFs[i].es[j];
			mm.FFs[i].es[0] = e0;
			mm.FFs[i].es[1] = e1;
			mm.FFs[i].es[2] = e2;
			mm.FFs[i].es[3] = e3;

			for (Int j = 0; j < 4; j++) mm.FEs[mm.FFs[i].es[j]].n_ffs.push_back(i);

			mm.FFs[i].vs.resize(4);
			mm.FFs[i].vs[0] = v0;
			mm.FFs[i].vs[1] = v1;
			mm.FFs[i].vs[2] = v2;
			mm.FFs[i].vs[3] = v3;

			for (Int j = 0; j < 4; j++) mm.FVs[mm.FFs[i].vs[j]].n_ffs.push_back(i);
		}

		return true;
	}

	void cuboid_extraction(Singularity &si, MetaMesh &mm, Mesh &mesh)
	{
		mm.FHs.clear();

		Int INVALID_F = (Int)mesh.Fs.size(), INVALID_H = (Int)mesh.Cs.size();
		std::vector<Int> f_flag(mesh.Fs.size(), INVALID_F), h_flag(mesh.Cs.size(), INVALID_H);
		std::vector<bool> ff_flag(mm.FFs.size(), false);

		for (Int i = 0; i < mm.FFs.size(); i++) for (Int j = 0; j < mm.FFs[i].ffs_grid.size(); j++)
			f_flag[mm.FFs[i].ffs_grid[j]] = i;

		while (true) {
			MetaMesh_H fh; fh.id = (Int)mm.FHs.size(); fh.Color_ID = -1;
			Int start_h = INVALID_H;
			for (Int i = 0; i < h_flag.size(); i++) if (h_flag[i] == INVALID_H) { start_h = i; break; }
			if (start_h == INVALID_H) break;

			std::queue<Int> h_pool; h_pool.push(start_h);
			while (!h_pool.empty()) {
				start_h = h_pool.front(); h_pool.pop();
				if (h_flag[start_h] != INVALID_H) continue;
				h_flag[start_h] = fh.id;
				fh.hs_grid.push_back(start_h);
				for (Int i = 0; i < 6; i++) {
					Int fid = mesh.Cs[start_h].fs[i];
					if (f_flag[fid] != INVALID_F) {
						if (!ff_flag[f_flag[fid]]) {
							ff_flag[f_flag[fid]] = true;
							fh.fs.push_back(f_flag[fid]);
						}
						continue;
					}
					for (Int j = 0; j < mesh.Fs[fid].n_cs.size(); j++) {
						Int hid = mesh.Fs[fid].n_cs[j];
						if (h_flag[hid] != INVALID_H) continue;
						h_pool.push(hid);
					}
				}
			}
			mm.FHs.push_back(fh);
			for (Int k = 0; k < fh.fs.size(); k++) ff_flag[fh.fs[k]] = false;
		}

		for (Int i = 0; i < mm.FHs.size(); i++) {
			mm.FHs[i].es.reserve(12);
			for (Int j = 0; j < mm.FHs[i].fs.size(); j++) {
				Int fid = mm.FHs[i].fs[j];
				for (Int k = 0; k < 4; k++) mm.FHs[i].es.push_back(mm.FFs[fid].es[k]);
			}
			std::sort(mm.FHs[i].es.begin(), mm.FHs[i].es.end());
			mm.FHs[i].es.erase(std::unique(mm.FHs[i].es.begin(), mm.FHs[i].es.end()), mm.FHs[i].es.end());
			mm.FHs[i].vs = mm.FFs[mm.FHs[i].fs[0]].vs;
			std::vector<Int> vs = mm.FFs[mm.FHs[i].fs[0]].vs;
			std::sort(vs.begin(), vs.end());
			short cors_f = 1;
			for (Int j = 1; j < 6; j++) {
				std::vector<Int> vsj = mm.FFs[mm.FHs[i].fs[j]].vs;
				std::sort(vsj.begin(), vsj.end());
				std::vector<Int> common_vs;
				std::set_intersection(vs.begin(), vs.end(), vsj.begin(), vsj.end(), std::back_inserter(common_vs));
				if (common_vs.size())continue; else { cors_f = j; break; }
			}
			vs = mm.FFs[mm.FHs[i].fs[cors_f]].vs;
			for (Int j = 0; j < 4; j++) {
				std::vector<Int> nvs = mm.FVs[mm.FHs[i].vs[j]].n_fvs;
				for (Int k = 0; k < nvs.size(); k++)
					if (std::find(vs.begin(), vs.end(), nvs[k]) != vs.end()) {
						mm.FHs[i].vs.push_back(nvs[k]); break;
					}
			}

			for (Int j = 0; j < 8; j++)mm.FVs[mm.FHs[i].vs[j]].n_fhs.push_back(i);
			for (Int j = 0; j < 12; j++)mm.FEs[mm.FHs[i].es[j]].n_fhs.push_back(i);
			for (Int j = 0; j < 6; j++)mm.FFs[mm.FHs[i].fs[j]].n_fhs.push_back(i);
		}

		for (auto &v : mm.FVs) v.boundary = false;
		for (auto &e : mm.FEs) e.boundary = false;
		for (auto &f : mm.FFs) {
			f.boundary = false;
			if (f.n_fhs.size() == 1) {
				f.boundary = true;
				for (auto vid : f.vs) mm.FVs[vid].boundary = true;
				for (auto eid : f.es) mm.FEs[eid].boundary = true;
			}
		}
	}

	void singularity_meta_mesh(Singularity &si, MetaMesh &mm, Mesh &mesh) {

		for (auto &fv : mm.FVs) {
			if (mesh.Vs[fv.hid].svid == -1) fv.svid = -1;
			else fv.svid = mesh.Vs[fv.hid].svid;
		}
		for (auto &fe : mm.FEs) {
			if ((!fe.boundary && fe.n_fhs.size() == Interior_RegularE) ||
				(fe.boundary && fe.n_fhs.size() == Boundary_RegularE))
				fe.singular = false;
			else
				fe.singular = true;
		}
	}

	void coloring(MetaMesh &frame)
	{
		int Color_id = 0; int Color_Total = 10;
		for (Int i = 0; i < frame.FHs.size(); i++) {
			if (frame.FHs[i].Color_ID == -1) {
				vector<int> ids;
				for (auto fid : frame.FHs[i].fs) {
					for (auto hid : frame.FFs[fid].n_fhs) {
						if (hid != i && frame.FHs[hid].Color_ID != -1)
						{
							ids.push_back(frame.FHs[hid].Color_ID);
						}
					}
				}
				while (true) {
					if (find(ids.begin(), ids.end(), Color_id) != ids.end()) {
						Color_id++;
						Color_id = Color_id % Color_Total;
					}
					else
						break;
				}
				frame.FHs[i].Color_ID = Color_id;
			}
		}
	}

	bool extract_metamesh(Mesh &mesh, Singularity &si, MetaMesh &mm)
	{
		mm.FVs.clear(); mm.FEs.clear(); mm.FFs.clear(); mm.FHs.clear();

		node_edge_extraction(si, mm, mesh);
		if (!face_extraction(si, mm, mesh)) return false;
		cuboid_extraction(si, mm, mesh);

		singularity_meta_mesh(si, mm, mesh);
		coloring(mm);

		return true;
	}

	bool extract_sheet(const Mesh &m, const vector<Int> &cand, vector<vector<Int>> &output, vector<std::pair<Int, Int>> &edge2sheet)
	{
		//only works for polycube-based hex mesh
		Int INVALID_ID = (Int)-1;

		vector<int> face_color(m.Fs.size(), -1), cand_flag(m.Fs.size(), -1);
		
		for (size_t i = 0; i < cand.size(); i++)
		{
			cand_flag[cand[i]] = 1;
		}

		output.clear();
		edge2sheet.clear();
		edge2sheet.resize(m.Es.size(), std::pair<Int, Int>(INVALID_ID, INVALID_ID));
		set<Int> cand_set(cand.begin(), cand.end());
		while (!cand_set.empty())
		{
			std::cout << "candidate set size: " << cand_set.size() << std::endl;
			//find boundary id;
			vector<Int> one_sheet;
			Int first_id = *cand_set.begin();
			face_color[first_id] = 1;
			queue<Int> q;
			q.push(first_id);
			while (!q.empty())
			{
				Int top = q.front(); q.pop();
				//iterate over 4 edges
				one_sheet.push_back(top);
				for (size_t i = 0; i < 4; i++)
				{
					Int edgeid = m.Fs[top].es[i];
					Int edgedegree = (Int)m.Es[edgeid].n_cs.size();

					//edge belong to sheet
					if (edge2sheet[edgeid].first == INVALID_ID)
						edge2sheet[edgeid].first = (Int)output.size();
					else if (edge2sheet[edgeid].first != output.size())
					{
						assert(edge2sheet[edgeid].second == INVALID_ID || edge2sheet[edgeid].second == output.size());
						edge2sheet[edgeid].second = (Int)output.size();
					}

					Int select_face = INVALID_ID;
					switch (edgedegree)
					{
					case 1: //edge turn here, not consider anymore
						break;
					case 2:
						if (m.Fs[top].boundary)
						{
							//find another boundary
							vector<int> twoboundary;
							for (size_t j = 0; j < 3; j++)
							{
								if (m.Fs[m.Es[edgeid].n_fs[j]].boundary)
									twoboundary.push_back(m.Es[edgeid].n_fs[j]);
							}
							assert(twoboundary.size() == 2);
							if (top == twoboundary[0] && cand_flag[twoboundary[1]] == 1)
								//q.push(twoboundary[1]);
								select_face = twoboundary[1];
							else if (top == twoboundary[1] && cand_flag[twoboundary[0]] == 1)
								select_face = twoboundary[0];
								//q.push(twoboundary[0]);

						}

						break;
					case 3:
						if (m.Fs[top].boundary)
						{
							//select a face which is not the boundary
							const vector<Int> &n_fs = m.Es[edgeid].n_fs;
							for (size_t j = 0; j < n_fs.size(); j++)
							{
								bool face1_flag = true;
								Int face_id = n_fs[j];
								if (m.Fs[face_id].boundary || cand_flag[face_id] == -1) continue;
								for (size_t k = 0; k < m.Fs[face_id].n_cs.size(); k++)
								{
									Int cell_id = m.Fs[face_id].n_cs[k];
									const vector<Int> &cell_faces = m.Cs[cell_id].fs;
									for (size_t p = 0; p < cell_faces.size(); p++)
									{
										if (cell_faces[p] == top)
										{
											face1_flag = false;
											break;
										}
									}
									if (!face1_flag)
										break;
								}

								if (face1_flag)
								{
									//q.push(n_fs[j]);
									select_face = n_fs[j];
									break;
								}
							}

							break;
						}
						//else go to default
						
					default:
						const vector<Int> &n_fs = m.Es[edgeid].n_fs;
						for (size_t j = 0; j < n_fs.size(); j++)
						{
							bool face1_flag = true;
							Int face_id = n_fs[j];
							if (face_id == top || cand_flag[face_id] == -1) continue;
							for (size_t k = 0; k < m.Fs[face_id].n_cs.size(); k++)
							{
								Int cell_id = m.Fs[face_id].n_cs[k];
								const vector<Int> &cell_faces = m.Cs[cell_id].fs;
								for (size_t p = 0; p < cell_faces.size(); p++)
								{
									if (cell_faces[p] == top)
									{
										face1_flag = false;
										break;
									}
								}
								if (!face1_flag)
									break;
							}
							
							if (face1_flag)
							{
								//q.push(n_fs[j]);
								select_face = n_fs[j];
								break;
							}
						}

						break;
					}

					if (select_face != INVALID_ID && face_color[select_face] != 1)
					{
						q.push(select_face);
						face_color[select_face] = 1;
					}
						
				}
			}
			output.push_back(one_sheet);
			//erase set
			for (size_t i = 0; i < one_sheet.size(); i++)
			{
				cand_set.erase(one_sheet[i]);
			}
		}

		return true;
	}

	Int face_cell_consistency(const Face &f, const Cell &c)
	{
		//return odd if face[bf_id] is towards outside, even otherwise
		//return value / 2 = index of the cell the face belongs to
		
		Int outside_flag = 1;
		Int face_idx = (Int)-1;

		const vector<Int> &allverts = c.vs;

		vector<vector<Int>> allfaces;
		for (size_t i = 0; i < 6; i++)
		{
			vector<Int> oneface;
			for (size_t j = 0; j < 4; j++)
			{
				oneface.push_back(allverts[hex_face_outside[i][j]]);
			}
			allfaces.push_back(oneface);
		}

		const vector<Int> bfv = f.vs;
		vector<Int> bfv_sort;
		for (size_t i = 0; i < bfv.size(); i++)
		{
			bfv_sort.push_back(bfv[i]);
		}
		sort(bfv_sort.begin(), bfv_sort.end());

		for (Int i = 0; i < 6; i++)
		{
			//find target face
			vector<Int> oneface_sort = allfaces[i];
			sort(oneface_sort.begin(), oneface_sort.end());
			Int equalcount = 0;
			for (size_t j = 0; j < 4; j++)
			{
				if (oneface_sort[j] == bfv_sort[j])
					equalcount++;
			}
			if (equalcount == 4)
			{
				//face found
				face_idx = i;
				for (size_t j = 0; j < 4; j++)
				{
					for (size_t k = 0; k < 4; k++)
					{
						if (bfv[j] == allfaces[i][k])
						{
							Int postj = (j + 1) % 4;
							Int postk = (k + 1) % 4;
							if (bfv[postj] == allfaces[i][postk])
							{
								outside_flag = 1;
								return (2 * face_idx + outside_flag);
								//return 1;
							}
							else
							{
								outside_flag = 0;
								return (2 * face_idx + outside_flag);
								//return 0;
							}
						}
					}
				}
			}

		}

		std::cout << "searching failed" << std::endl;
		
		return (2 * face_idx + outside_flag);
	}

	Int get_boundary_face_dir(const Mesh &m, Int bf_id)
	{
		//return 1 if face[bf_id] is towards outside, 0 otherwise
		assert(m.Fs[bf_id].boundary == true);
		Int cellid = m.Fs[bf_id].n_cs[0];
		return face_cell_consistency(m.Fs[bf_id], m.Cs[cellid]) % 2;
	}
	
	bool get_neighbor_face_consistency(const Mesh &m, Int id0, Int id1)
	{
		std::vector<Int> hf0_vert_idx = m.Fs[id0].vs, hf1_vert_idx = m.Fs[id1].vs; //no ref

		std::reverse(hf1_vert_idx.begin(), hf1_vert_idx.end());

		std::pair<Int, Int> begin_idx(-1, -1);
		for (Int i = 0; i < hf0_vert_idx.size(); i++)
		{

			for (Int j = 0; j < hf1_vert_idx.size(); j++)
			{
				if (hf0_vert_idx[i] == hf1_vert_idx[j])
				{
					Int post0 = (Int)(i + 1) % hf0_vert_idx.size();
					Int post1 = (Int)(j + 1) % hf1_vert_idx.size();
					Int pre0 = (Int)(i - 1 + hf0_vert_idx.size()) % hf0_vert_idx.size();
					Int pre1 = (Int)(j - 1 + hf1_vert_idx.size()) % hf1_vert_idx.size();
					if (hf0_vert_idx[post0] == hf1_vert_idx[post1] || hf0_vert_idx[pre0] == hf1_vert_idx[pre1])
						return true;
					else
						return false;
				}
			}
		}

		std::cout << "no common verts! ! !" << std::endl;

		return false;
	}

	void get_consistent_rotation(const Mesh &m, const vector<Int> &sheet, vector<Int>& rotationflag, bool refineflag)
	{
		//space of rotation_flag should be pre-assigned
		//adjust rotation flag so that each face is towards outside
		//1: default order 0: reverse order 
		Int INVALID_ID = (Int)-1;

		Int boundaryid = INVALID_ID;
		Int nf = (Int)rotationflag.size();
		nf = nf < m.Fs.size() ? (Int)m.Fs.size() : nf;
		vector<bool> sheetflag(nf, false);
		queue<Int> q;

		for (size_t i = 0; i < sheet.size(); i++)
		{
			sheetflag[sheet[i]] = true;
		}

		if (!refineflag)
		{
			//rotationflag.clear();
			
			//rotationflag.resize(nf, INVALID_ID);

			//decide each sheet
			
			for (size_t i = 0; i < sheet.size(); i++)
			{
				if (m.Fs[sheet[i]].boundary)
				{
					boundaryid = sheet[i];
					break;
				}
			}

			if (boundaryid != INVALID_ID)
			{
				//decide rotation for boundary id
				rotationflag[boundaryid] = get_boundary_face_dir(m, boundaryid);
			}
			else
			{
				//select the face that close to a boundary cell
				boundaryid = sheet[0];
				rotationflag[boundaryid] = 1;
			}
		}
		else
		{
			//just refine
			for (size_t i = 0; i < sheet.size(); i++)
			{
				if (rotationflag[sheet[i]] != INVALID_ID)
				{
					boundaryid = sheet[i];
					break;
				}
			}
		}
		//assert(boundaryid != INVALID_ID);
		//assertion not always correct
		if (boundaryid == INVALID_ID) boundaryid = sheet[0];
		q.push(boundaryid);
		vector<Int> facecolor(nf, 0);
		facecolor[boundaryid] = 1;

		while (!q.empty())
		{
			Int front = q.front(); q.pop();
			
			//find nb of front
			const vector<Int> &ne = m.Fs[front].es;

			for (size_t i = 0; i < ne.size(); i++)
			{
				const vector<Int> &nf = m.Es[ne[i]].n_fs;
				for (size_t j = 0; j < nf.size(); j++)
				{
					if (nf[j] == front || !sheetflag[nf[j]] || facecolor[nf[j]] == 1) continue;

					bool consistency = get_neighbor_face_consistency(m, front, nf[j]);

					if (consistency)
					{
						rotationflag[nf[j]] = rotationflag[front];
					}
					else
					{
						rotationflag[nf[j]] = 1 - rotationflag[front];
					}
					q.push(nf[j]);
					facecolor[nf[j]] = 1;
				}
			}

		}



	}

	void compute_cell_quality(Mesh &m, int type)
	{
		//compute hex quality for all data type£¬ store them in cell.quality
		//0: scaled jacobian
		//1: volume
		//2: quality

		size_t nc = m.Cs.size();
		for (size_t i = 0; i < nc; i++)
		{
			std::vector<CVec<double, 3>> elem_pts;
			for (size_t j = 0; j < 8; j++)
			{
				elem_pts.push_back(m.Vs[m.Cs[i].vs[j]].v);
			}
			
			Hex hex(elem_pts);
			m.Cs[i].quality = hex.compute_quality(type);
			//assert(m.Cs[i].quality > 0.0);
			
			/*if (i < 100)
				std::cout << m.Cs[i].quality << std::endl;*/
			
		}

	}

}