#include <cassert>

namespace ig
{
    namespace aabb_impl
    {
        //////////////////////////////////////////////////////////////
        template <typename T>
        AABB_Segment_Tree_Node<T>::AABB_Segment_Tree_Node(const std::vector<Vector3>& vertexList, const std::vector<ptrdiff_t>& segment_id_list, AABBTreeInfo& treeInfo, int depth)
        {
            left = nullptr;
            right = nullptr;
            m_nNumVerts = vertexList.size();
            BuildTree(vertexList, segment_id_list, treeInfo, depth);
        }
        //////////////////////////////////////////////////////////////
        template <typename T>
        AABB_Segment_Tree_Node<T>::~AABB_Segment_Tree_Node()
        {
            if (left) delete left;
            if (right) delete right;
        }
        //////////////////////////////////////////////////////////////
        template <typename T>
        int AABB_Segment_Tree_Node<T>::FindBestAxis(const std::vector<Vector3>& vertexList)
        {
            size_t vertexNum = vertexList.size();

            //divide this box into two boxes - pick a better axis
            int iAxis = 0;
            int iAxisResult[3]; //stores how close end result is, the lower the better

            Vector3 center = m_box.getCenter();

            for (iAxis = 0; iAxis < 3; iAxis++)
            {
                int ileft = 0, iright = 0;
                Vector3 v[2];
                int count = 0;
                for (int i = 0; i < vertexNum; i++)
                {
                    v[count] = vertexList[i];
                    if (count == 1)
                    {
                        T MidPoint[3];
                        MidPoint[0] = (v[0][0] + v[1][0]) / (T)2.0;
                        MidPoint[1] = (v[0][1] + v[1][1]) / (T)2.0;
                        MidPoint[2] = (v[0][2] + v[1][2]) / (T)2.0;

                        if (MidPoint[iAxis] <= center[iAxis])
                        {
                            ileft++;
                        }
                        else
                        {
                            iright++;
                        }

                        count = 0;
                    }
                    else
                        count++;
                } // vertices

                iAxisResult[iAxis] = abs(ileft - iright);
            } //axis

            int index = 0;
            int result = iAxisResult[0];
            for (int i = 1; i < 3; i++)
            {
                if (iAxisResult[i] < result)
                {
                    result = iAxisResult[i];
                    index = i;
                }
            }

            return index;
        }
        //////////////////////////////////////////////////////////////
        template <typename T>
        void AABB_Segment_Tree_Node<T>::SetBounds(const aabb_impl::AABB<T>& box)
        {
            min = box.min;
            max = box.max;
            m_box = box;
        }
        //////////////////////////////////////////////////////////////
        template <typename T>
        void AABB_Segment_Tree_Node<T>::BuildTree(const std::vector<Vector3>& vertexList, const std::vector<ptrdiff_t>& segment_id_list, AABBTreeInfo& treeInfo, int depth)
        {
            assert(segment_id_list.size() * 2 == vertexList.size());
            ptrdiff_t vertexNum = vertexList.size();

            // Build the node bounding box based from the triangle list
            aabb_impl::AABB<T> Box(&vertexList[0], vertexNum);

            //debug box bounds
            SetBounds(Box);

            if (depth + 1 > treeInfo.curr_max_depth)
                treeInfo.curr_max_depth = depth + 1;

            bool bMakeChildren = false;

            if (vertexNum > treeInfo.min_vertices && depth < treeInfo.max_tree_depth)
            {
                bMakeChildren = true;
            }

            if (bMakeChildren)
            {
                // Find the longest axii of the node's box
                int iAxis = FindBestAxis(vertexList);
                //Log("Longest axis: %d\n", iAxis);

                //Get the Arrays for min, max dimensions
                //float *min = &Box.min[0];
                //float *max = &Box.max[0];

                //get center of the box for longest axis
                Vector3 center = Box.getCenter();

                int count = 0;
                Vector3 v[2];
                std::vector<Vector3> leftSide;
                std::vector<Vector3> rightSide;
                std::vector<ptrdiff_t> leftSide_segment_id;
                std::vector<ptrdiff_t> rightSide_segment_id;

                int leftCount = 0, rightCount = 0; //debug

                //btw, things that could go wrong- if the mesh doesn't
                //send over the triangles correctly, then you might see
                //huge boxes that are misaligned (bad leaves).
                //things to check is making sure the vertex buffer is
                //correctly aligned along the adjancey buffers, etc
                for (int i = 0; i < vertexNum; i++)
                {
                    v[count] = vertexList[i];

                    if (count == 1)
                    {
                        T MidPoint[3];
                        MidPoint[0] = (v[0][0] + v[1][0]) / (T)2.0;
                        MidPoint[1] = (v[0][1] + v[1][1]) / (T)2.0;
                        MidPoint[2] = (v[0][2] + v[1][2]) / (T)2.0;
                        //WVector faceCenter(x,y,z);

                        if (MidPoint[iAxis] <= center[iAxis]) //fSplit
                        {
                            //Store the verts to the left.
                            leftSide.push_back(v[0]);
                            leftSide.push_back(v[1]);
                            leftSide_segment_id.push_back(segment_id_list[i / 2]);
                            leftCount++;
                        }
                        else
                        {
                            //Store the verts to the right.
                            rightSide.push_back(v[0]);
                            rightSide.push_back(v[1]);
                            rightSide_segment_id.push_back(segment_id_list[i / 2]);
                            rightCount++;
                        }

                        count = 0;
                    }
                    else
                        count++;
                }

                if (treeInfo.m_bPrune && (leftCount == 0 || rightCount == 0))
                {
                    //okay, now it's time to cheat. we couldn't use
                    //the best axis to split the
                    //box so now we'll resort to brute force hacks....
                    leftSide.clear();
                    rightSide.clear();
                    leftSide_segment_id.clear();
                    rightSide_segment_id.clear();

                    ptrdiff_t leftMaxIndex = vertexNum / 2; //left side

                    count = 0;
                    for (int i = 0; i < vertexNum; i++)
                    {
                        v[count] = vertexList[i];
                        if (count == 1)
                        {
                            if (i < leftMaxIndex)
                            {
                                //left node
                                leftSide.push_back(v[0]);
                                leftSide.push_back(v[1]);
                                leftSide_segment_id.push_back(segment_id_list[i / 2]);
                            }
                            else
                            {
                                rightSide.push_back(v[0]);
                                rightSide.push_back(v[1]);
                                rightSide_segment_id.push_back(segment_id_list[i / 2]);
                            }

                            count = 0;
                        }
                        else
                            count++;
                    }
                }

                if (leftSide.size() > 0 && rightSide.size() > 0)
                {
                    assert(leftSide.size() % 2 == 0);
                    assert(rightSide.size() % 2 == 0);

                    //Build child nodes
                    if (leftSide.size() > 0)
                    {
                        treeInfo.left_children++;
                        left = new AABB_Segment_Tree_Node(leftSide, leftSide_segment_id, treeInfo, depth + 1);
                    }
                    if (rightSide.size() > 0)
                    {
                        treeInfo.right_children++;
                        right = new AABB_Segment_Tree_Node(rightSide, rightSide_segment_id, treeInfo, depth + 1);
                    }
                }
                else
                {
                    //should never happen
                    bMakeChildren = false;
                }
            }

            if (!bMakeChildren)
            {
                //Store the data directly if you want....
                for (int i = 0; i < vertexList.size(); i++)
                {
                    m_pVerts.push_back(vertexList[i]);
                }
                for (int i = 0; i < segment_id_list.size(); i++)
                {
                    segment_id.push_back(segment_id_list[i]);
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////
    template <typename T>
    AABB_Segment_Tree<T>::AABB_Segment_Tree(const std::vector<Vector3>& vertices_list)
    {
        treeInfo.max_tree_depth = 100;
        treeInfo.min_vertices = 50;
        treeInfo.curr_max_depth = 0;
        treeInfo.m_bPrune = true;
        std::vector<ptrdiff_t> segment_id_list;
        segment_id_list.reserve(vertices_list.size() / 2);
        for (int i = 0; i < vertices_list.size() / 2; i++)
            segment_id_list.push_back(i);
        vertices_list_copy.assign(vertices_list.begin(), vertices_list.end());
        root = new aabb_impl::AABB_Segment_Tree_Node<T>(vertices_list, segment_id_list, treeInfo, 0);
        kdTree = 0;
        if (Enable_AABB_Segment_KDTree)
        {
            std::vector<Vector3> v_list;
            std::vector<ptrdiff_t> local_segment_id_list;
            //sample some points
            for (int i = 0; i < vertices_list.size() / 2; i++)
            {
                v_list.push_back(vertices_list[2 * i + rand() % 2]);
				local_segment_id_list.push_back(i);
                v_list.push_back((vertices_list[2 * i] + vertices_list[2 * i + 1]) / 2.0);
				local_segment_id_list.push_back(i);
            }
            kdTree = new KD_Tree<T, 3>(v_list, local_segment_id_list);
        }
    }
    //////////////////////////////////////////////////////////////
    template <typename T>
    AABB_Segment_Tree<T>::~AABB_Segment_Tree()
    {
        if (root) delete root;
        if(kdTree) delete kdTree;
    }
    //////////////////////////////////////////////////////////////
    template <typename T>
    T AABB_Segment_Tree<T>::point_AABB_distance(const Vector3& p, const aabb_impl::AABB<T>& box)
    {
        T minx, miny, minz, maxx, maxy, maxz;
        T x, y, z;
        T mindist;
        x = p[0];
        y = p[1];
        z = p[2];
        minx = box.min[0];
        miny = box.min[1];
        minz = box.min[2];
        maxx = box.max[0];
        maxy = box.max[1];
        maxz = box.max[2];
        Vector3 corner[8];
        if (x >= minx && x <= maxx && y >= miny && y <= maxy && z >= minz && z <= maxz)
        {
            //inside AABB
            return -1;
        }
        mindist = 0;
        if (x < minx)
        {
            mindist += (x - minx)*(x - minx);
        }
        else if (x>maxx)
        {
            mindist += (x - maxx)*(x - maxx);
        }

        if (y < miny)
        {
            mindist += (miny - y)*(miny - y);
        }
        else if (y>maxy)
        {
            mindist += (y - maxy)*(y - maxy);
        }

        if (z < minz)
        {
            mindist += (minz - z)*(minz - z);
        }
        else if (z>maxz)
        {
            mindist += (z - maxz)*(z - maxz);
        }

        return sqrt(mindist);
    }
    //////////////////////////////////////////////////////////////
    template <typename T>
    T AABB_Segment_Tree<T>::point_line_distance(const Vector3& p, const Vector3& a, const Vector3& b, Vector3& nearestP, T ref_dist)
    {
        Vector3 ab = b - a;
        Vector3 pa = a - p;
        Vector3 p_inline;
        T point_line_dist;
        T rst_dist;
        ab.Normalize();
        point_line_dist = (pa - pa.Dot(ab)*ab).Length();
        p_inline = p + pa - pa.Dot(ab)*ab;

        if (ref_dist >= 0 && point_line_dist > ref_dist)
        {
            //reduction of computation
            return point_line_dist;
        }

        if (inside_segment(p_inline, a, b))
        {
            nearestP = p_inline;
            rst_dist = point_line_dist;
        }
        else
        {
            nearestP = a;
            rst_dist = (a - p).Length();
            if ((b - p).Length() < rst_dist)
            {
                rst_dist = (b - p).Length();
                nearestP = b;
            }
        }
        return rst_dist;
    }
    //////////////////////////////////////////////////////////////
    template <typename T>
    bool AABB_Segment_Tree<T>::inside_segment(const Vector3& p, const Vector3& a, const Vector3& b)
    {
        Vector3 pa, pb;
        pa = a - p;
        pb = b - p;
        pa.Normalize();
        pb.Normalize();
        if (fabs(fabs(pa.Dot(pb)) - 1) > (T)1e-6)
            return false;
        if (pa.Dot(pb) < -(T)0.5)
            return true;
        return false;
    }
    //////////////////////////////////////////////////////////////
    template <typename T>
    T AABB_Segment_Tree<T>::search_entire_node(aabb_impl::AABB_Segment_Tree_Node<T>* node, const Vector3& p, Vector3& globalVector3, T& global_min_dist, ptrdiff_t& global_Segment_id)
    {
        T min_dist, temp_rst;
        ptrdiff_t cur_segment_id;
        Vector3 min_point;
        Vector3 rst;
        if (node == nullptr)
            return -1;
        if (node->m_nNumVerts < 3)
            return -1;

        min_dist = point_line_distance(p, node->m_pVerts[0], node->m_pVerts[1], min_point, -1);
        cur_segment_id = node->segment_id[0];
        for (int i = 1; i < node->m_nNumVerts / 2; i++)
        {
            temp_rst = point_line_distance(p, node->m_pVerts[2 * i], node->m_pVerts[2 * i + 1], rst, min_dist);
            if (temp_rst < min_dist)
            {
                min_dist = temp_rst;
                min_point = rst;
                cur_segment_id = node->segment_id[i];
            }
        }
        if (min_dist < global_min_dist)
        {
            global_min_dist = min_dist;
            globalVector3 = min_point;
            global_Segment_id = cur_segment_id;
        }
        return global_min_dist;
    }
    //////////////////////////////////////////////////////////////
    template <typename T>
    void AABB_Segment_Tree<T>::find_Nearest_Point(const Vector3& p, Vector3& nearestP, ptrdiff_t* segment_id)
    {
        T global_min_dist;
        Vector3 globalVector3;
        ptrdiff_t global_Segment_id;
        if (Enable_AABB_Segment_KDTree)
        {
            ptrdiff_t cur_segment_id;
            kdTree->find_Nearest_Point(p, globalVector3, &cur_segment_id);
            global_Segment_id = cur_segment_id;
            global_min_dist = (globalVector3 - p).Length();
        }
        else
        {
            globalVector3 = vertices_list_copy[0];
            global_min_dist = (vertices_list_copy[0] - p).Length();
            global_Segment_id = 0;
        }

        Traverse_Search(root, p, globalVector3, global_min_dist, global_Segment_id);
        nearestP = globalVector3;
        if (segment_id != nullptr)
        {
            if (global_Segment_id == -1)
            {
                //should never goes here
                for (int i = 0; i < vertices_list_copy.size(); i++)
                    if ((vertices_list_copy[i] - globalVector3).Length() < 1e-10)
                    {
                        global_Segment_id = i / 2;
                        break;
                    }
            }
            *segment_id = global_Segment_id;
        }
       // return global_min_dist;
    }
    //////////////////////////////////////////////////////////////
    template <typename T>
    void AABB_Segment_Tree<T>::Traverse_Search(aabb_impl::AABB_Segment_Tree_Node<T>* node, const Vector3& p, Vector3& globalVector3, T& global_min_dist, ptrdiff_t& global_Segment_id)
    {
        if (node == nullptr)
            return;
        if (point_AABB_distance(p, node->m_box) > global_min_dist)
            return;
        if (node->left == nullptr && node->right == nullptr)
        {
            //search entire node
            search_entire_node(node, p, globalVector3, global_min_dist, global_Segment_id);
            return;
        }
        if (node->left != nullptr)
            Traverse_Search(node->left, p, globalVector3, global_min_dist, global_Segment_id);
        if (node->right != nullptr)
            Traverse_Search(node->right, p, globalVector3, global_min_dist, global_Segment_id);
    }
    //////////////////////////////////////////////////////////////

}