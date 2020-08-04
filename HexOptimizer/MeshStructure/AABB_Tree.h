#pragma once
/*
This Two Class AABB_Tree and AABB_Segment_Tree is for fast point-(triangle soup / line segment soup) nearest distance query
There is only one constructor in each class.

AABB_Tree(std::vector<Vertex> &vertices_list);
AABB_Segment_Tree(std::vector<Vertex> &vertices_list);

In AABB_Tree vertices_list is a vector size=3n store n triangle vertices(tri_1_p1,tri_1_p2,tri_1_p3,tri_2_p1...tri_n_p3)
In AABB_Segment_Tree vertices_list is a vector size=2n store n segment endpoints(segment_1_p1,segment_1_p2,segment_2_p1...segment_n_p2)

And there is only one public function for user to use;

float findNearstPoint(Vector3 p, Vector3 &nearestP,int *face_id=nullptr);
float findNearstPoint(Vector3 p, Vector3 &nearestP, int *segment_id = nullptr);

where p is query point, and nearest point will be stored in nearestP, if face_id/segment_id not equal to nullptr, the corresponding primitive id of nearest point will be
stored in face_id/segment_id(id order follow the vertices_list storage order form 0 to n-1)
*/

#include "KD_Tree.h"

namespace ig
{
#define Enable_AABB_KDTree 1
#define Enable_AABB_Segment_KDTree 1
    /*
    	There is an option for user to decide whether to use KD_TREE or not.
    	KD_Tree will help to accelerate query with some extra preprocessing time and extra memory to store KD_TREE.
    */

    namespace aabb_impl
    {
        template <typename T>
        class AABB
        {
            typedef TinyVector<T, 3> Vector3;
        public:
            AABB();
            AABB(const Vector3* vertex_list, ptrdiff_t num);
            AABB(const std::vector<Vector3>& vertices_list);
            Vector3 getCenter();
        public:
            Vector3 min, max;
        };

        struct AABBTreeInfo
        {
            //heuristics
            int max_tree_depth;		//max depth the tree can reach
            //min number of vertices each leaf can store
            int min_vertices;

            int curr_max_depth;
            //ensures that an AABB is never generated that
            //is over min_vertices. Normally, this algorithm
            //is not required because the best axis
            //algorithm normally produces a perfectly balanced tree.
            bool m_bPrune;
            int left_children, right_children;
        };

        template <typename T>
        class AABBNode
        {
            typedef TinyVector<T, 3> Vector3;
        public:
            AABBNode(const std::vector<Vector3>& vertexList, const std::vector<ptrdiff_t>& face_id_list, AABBTreeInfo& treeInfo, int depth);
            ~AABBNode();
        protected:
            int FindBestAxis(const std::vector<Vector3>& vertexList);

            void SetBounds(const AABB<T>& box);

            void BuildTree(const std::vector<Vector3>& vertexList, const std::vector<ptrdiff_t>& face_id_list, AABBTreeInfo& treeInfo, int depth);

        public:
            Vector3 min, max;
            AABBNode<T>* left, *right;
            std::vector<Vector3> m_pVerts;	//vertices stored within this node
            std::vector<ptrdiff_t> face_id;
            ptrdiff_t m_nNumVerts; //num of vertices stored (3->n)
            AABB<T> m_box;
        };
    }


    template <typename T>
    class AABB_Tree
    {
    public:
		typedef T Real;
        typedef TinyVector<T, 3> Vector3;
        AABB_Tree(const std::vector<Vector3>& vertices_list);
        ~AABB_Tree();
        void find_Nearest_Point(const Vector3& p, Vector3& nearestP, ptrdiff_t* face_ID = nullptr);
        void find_Nearest_Point_by_traverse(const Vector3& p, Vector3& nearestP, ptrdiff_t* face_ID = nullptr);
		Vector3 get_normal(ptrdiff_t face_ID);
    protected:
        T point_line_distance(const Vector3& p, const Vector3& a, const Vector3& b, Vector3& nearestP, T ref_dist);
        bool inside_segment(const Vector3& p, const Vector3& a, const Vector3& b);
        void Traverse(aabb_impl::AABBNode<T>* node);
        T point_AABB_distance(const Vector3& p, const aabb_impl::AABB<T>& box);
        void Traverse_Search(aabb_impl::AABBNode<T>* node, const Vector3& p, Vector3& globalVector3, T& global_min_dist, ptrdiff_t& globalFaceID);
		T search_entire_node(aabb_impl::AABBNode<T>* node, const Vector3& p, Vector3& globalVector3, T& global_min_dist, ptrdiff_t& globalFaceID);
        T point_tri_distance_refine(const Vector3& p, const Vector3& a, const Vector3& b, const Vector3& c, Vector3& nearestP, T ref_dist);
    private:
        KD_Tree<T, 3>* kdTree;
        aabb_impl::AABBNode<T>* root;
        aabb_impl::AABBTreeInfo treeInfo;
        ptrdiff_t leafNodeNum;
        ptrdiff_t all_face_num;
        std::vector<Vector3> vertices_list_copy;
    };

}
#include "AABB_Tree.inl"