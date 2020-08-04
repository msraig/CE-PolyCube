#include "TetStructure.h"
#include <set>
#include <map>
#include <vector>

namespace ig
{
	template <typename Real>
	class TetOperation
	{
	public:
		TetOperation() { tetstructure = NULL; };
		~TetOperation()
		{
			if (tetstructure)
			{
				//delete tetstructure;
				tetstructure = NULL;
			}
		}
		void set_tet_structure(TetStructure<Real>* tet)
		{
			if (tetstructure != tet)
			{
				delete tetstructure;
			}
			tetstructure = tet;
			
			//cell_matching construction
			cell_matching.clear();
			cell_matching.resize(tetstructure->tetras.size() * 4, 0);
		}

		bool split_edge(Tetrahedron<Real>* tet, unsigned int vi, unsigned int vj, bool split_boundary = true);

		bool split_inner_edge_two_boundary_point(TetStructure<Real> &tet_mesh);
		
		bool split_cell_four_boundary_point(TetStructure<Real> &tet_mesh);

		void get_vert_cell(TetStructure<Real>* tet_mesh, std::vector<CVec<Real, 3>> &points, std::vector<unsigned int> &indices);

		//if ori_flag equals to true, extract all verts
		void get_boundary_vert_face(TetStructure<Real>* tet_mesh, std::vector<Real> &points, std::vector<int> &faces, std::vector<int> &s2v, bool ori_vert_flag = false, std::map<std::vector<int>, int> *face2cellidx = NULL);

		void get_boundary_edge_set(const TetStructure<Real> &tet_mesh, std::set<std::pair<int, int>> &boundary_edge_set);

		void get_edge_local_id(const std::vector<Tetrahedron<Real>* > &tetras, std::pair<int, int> edge, int &tet_id, int localid[]);

	private:
		TetStructure<Real> *tetstructure;
		std::vector<int> cell_matching;
	};
}

