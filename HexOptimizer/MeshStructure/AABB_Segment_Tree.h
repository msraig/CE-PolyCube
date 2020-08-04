#pragma once
#include "AABB_Tree.h"

/*
 Notes in AABB_Tree.h
*/

namespace ig
{
    namespace aabb_impl
    {
        template <typename T>
        class AABB_Segment_Tree_Node
        {
            typedef TinyVector<T, 3> Vector3;
        public:
            AABB_Segment_Tree_Node(const std::vector<Vector3>& vertexList, const std::vector<ptrdiff_t>& segment_id_list, AABBTreeInfo& treeInfo, int depth);
            ~AABB_Segment_Tree_Node();
        private:
            int FindBestAxis(const std::vector<Vector3>& vertexList);
            void SetBounds(const aabb_impl::AABB<T>& box);
            void BuildTree(const std::vector<Vector3>& vertexList, const std::vector<ptrdiff_t>& segment_id_list, AABBTreeInfo& treeInfo, int depth);
        public:
            Vector3 min, max;
            AABB_Segment_Tree_Node<T>* left, *right;
            std::vector<Vector3> m_pVerts;	//vertices stored within this node
            std::vector<ptrdiff_t> segment_id;
            ptrdiff_t m_nNumVerts; //num of vertices stored (3->n)
            aabb_impl::AABB<T> m_box;
        };
    }

    template <typename T>
    class AABB_Segment_Tree
    {
    public:
		typedef T Real;
        typedef TinyVector<T, 3> Vector3;
        AABB_Segment_Tree(const std::vector<Vector3>& vertices_list);
        ~AABB_Segment_Tree();
        void find_Nearest_Point(const Vector3& p, Vector3& nearestP, ptrdiff_t* segment_id = nullptr);
    private:
        T point_AABB_distance(const Vector3& p, const aabb_impl::AABB<T>& box);
        T point_line_distance(const Vector3& p, const Vector3& a, const Vector3& b, Vector3& nearestP, T ref_dist);
        bool inside_segment(const Vector3& p, const Vector3& a, const Vector3& b);
        T search_entire_node(aabb_impl::AABB_Segment_Tree_Node<T>* node, const Vector3& p, Vector3& globalVector3, T& global_min_dist, ptrdiff_t& global_Segment_id);
        void Traverse_Search(aabb_impl::AABB_Segment_Tree_Node<T>* node, const Vector3& p, Vector3& globalVector3, T& global_min_dist, ptrdiff_t& global_Segment_id);

    private:
        aabb_impl::AABB_Segment_Tree_Node<T>* root;
        aabb_impl::AABBTreeInfo treeInfo;
        KD_Tree<T, 3>* kdTree;
        std::vector<Vector3> vertices_list_copy;
    };

}

#include "AABB_Segment_Tree.inl"