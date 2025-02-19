#include <cassert>

namespace ig
{
    namespace aabb_impl
    {
        //////////////////////////////////////////////////////////////
        template <typename T>
        AABB<T>::AABB()
        {
            min = Vector(0, 0, 0);
            max = Vector(0, 0, 0);
        }
        //////////////////////////////////////////////////////////////
        template <typename T>
        AABB<T>::AABB(const Vector *vertex_list, ptrdiff_t num)
        {
            min = max = vertex_list[0];

            for (ptrdiff_t i = 1; i < num; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    min[j] = std::min(min[j], vertex_list[i][j]);
                    max[j] = std::max(max[j], vertex_list[i][j]);
                }
            }
        }
        //////////////////////////////////////////////////////////////
        template <typename T>
        AABB<T>::AABB(const std::vector<Vector> &vertices_list)
        {
            min = max = vertices_list[0];

            for (ptrdiff_t i = 1; i < (ptrdiff_t)vertices_list.size(); i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    min[j] = std::min(min[j], vertices_list[i][j]);
                    max[j] = std::max(max[j], vertices_list[i][j]);
                }
            }
        }
        //////////////////////////////////////////////////////////////
        template <typename T>
        TinyVector<T, 3> AABB<T>::getCenter()
        {
            return (T)0.5 * (min + max);
        }
        //////////////////////////////////////////////////////////////
        template <typename T>
        AABBNode<T>::AABBNode(const std::vector<Vector> &vertexList, const std::vector<ptrdiff_t> &face_id_list, AABBTreeInfo &treeInfo, int depth)
        {
            left = nullptr;
            right = nullptr;
            m_nNumVerts = vertexList.size();
            BuildTree(vertexList, face_id_list, treeInfo, depth);
        }
        //////////////////////////////////////////////////////////////
        template <typename T>
        AABBNode<T>::~AABBNode()
        {
            if (left)
                delete left;
            if (right)
                delete right;
        }
        //////////////////////////////////////////////////////////////
        template <typename T>
        int AABBNode<T>::FindBestAxis(const std::vector<Vector> &vertexList)
        {
            size_t vertexNum = vertexList.size();

            // divide this box into two boxes - pick a better axis
            int iAxis = 0;
            int iAxisResult[3]; // stores how close end result is, the lower the better

            Vector center = m_box.getCenter();

            for (iAxis = 0; iAxis < 3; iAxis++)
            {
                int ileft = 0, iright = 0;
                Vector v[3];
                int count = 0;
                for (int i = 0; i < vertexNum; i++)
                {
                    v[count] = vertexList[i];
                    if (count == 2)
                    {
                        Vector faceCenter = (v[0] + v[1] + v[2]) / (T)3;
                        if (faceCenter[iAxis] <= center[iAxis])
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
            } // axis

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
        void AABBNode<T>::SetBounds(const AABB<T> &box)
        {
            min = box.min;
            max = box.max;
            m_box = box;
        }
        //////////////////////////////////////////////////////////////
        template <typename T>
        void AABBNode<T>::BuildTree(const std::vector<Vector> &vertexList, const std::vector<ptrdiff_t> &face_id_list, AABBTreeInfo &treeInfo, int depth)
        {
            int vertexNum = (int)vertexList.size();

            // Build the node bounding box based from the triangle list
            AABB<T> Box(vertexList);

            // debug box bounds
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
                // Log("Longest axis: %d\n", iAxis);

                // Get the Arrays for min, max dimensions
                // T* m_min = &Box.min[0];
                // T* m_max = &Box.max[0];

                // get center of the box for longest axis
                Vector center = Box.getCenter();

                int count = 0;
                Vector v[3];
                std::vector<Vector> leftSide;
                std::vector<Vector> rightSide;
                std::vector<ptrdiff_t> leftSide_face_id;
                std::vector<ptrdiff_t> rightSide_face_id;

                int leftCount = 0, rightCount = 0; // debug

                // btw, things that could go wrong- if the mesh doesn't
                // send over the triangles correctly, then you might see
                // huge boxes that are misaligned (bad leaves).
                // things to check is making sure the vertex buffer is
                // correctly aligned along the adjancey buffers, etc
                for (int i = 0; i < vertexNum; i++)
                {
                    v[count] = vertexList[i];

                    if (count == 2)
                    {
                        Vector faceCenter = (v[0] + v[1] + v[2]) / (T)3;

                        if (faceCenter[iAxis] <= center[iAxis]) // fSplit
                        {
                            // Store the verts to the left.
                            leftSide.emplace_back(v[0]);
                            leftSide.emplace_back(v[1]);
                            leftSide.emplace_back(v[2]);
                            leftSide_face_id.emplace_back(face_id_list[i / 3]);
                            leftCount++;
                        }
                        else
                        {
                            // Store the verts to the right.
                            rightSide.emplace_back(v[0]);
                            rightSide.emplace_back(v[1]);
                            rightSide.emplace_back(v[2]);
                            rightSide_face_id.emplace_back(face_id_list[i / 3]);
                            rightCount++;
                        }

                        count = 0;
                    }
                    else
                        count++;
                }

                if (treeInfo.m_bPrune && (leftCount == 0 || rightCount == 0))
                {
                    // okay, now it's time to cheat. we couldn't use
                    // the best axis to split the
                    // box so now we'll resort to brute force hacks....
                    leftSide.clear();
                    rightSide.clear();
                    leftSide_face_id.clear();
                    rightSide_face_id.clear();

                    int leftMaxIndex = vertexNum / 2; // left side

                    count = 0;
                    for (int i = 0; i < vertexNum; i++)
                    {
                        v[count] = vertexList[i];
                        if (count == 2)
                        {
                            if (i < leftMaxIndex)
                            {
                                // left node
                                leftSide.emplace_back(v[0]);
                                leftSide.emplace_back(v[1]);
                                leftSide.emplace_back(v[2]);
                                leftSide_face_id.emplace_back(face_id_list[i / 3]);
                            }
                            else
                            {
                                rightSide.emplace_back(v[0]);
                                rightSide.emplace_back(v[1]);
                                rightSide.emplace_back(v[2]);
                                rightSide_face_id.emplace_back(face_id_list[i / 3]);
                            }

                            count = 0;
                        }
                        else
                            count++;
                    }
                }

                if (leftSide.size() > 0 && rightSide.size() > 0)
                {
                    assert(leftSide.size() % 3 == 0);
                    assert(rightSide.size() % 3 == 0);

                    // Build child nodes
                    if (leftSide.size() > 0)
                    {
                        treeInfo.left_children++;
                        left = new AABBNode(leftSide, leftSide_face_id, treeInfo, depth + 1);
                    }
                    if (rightSide.size() > 0)
                    {
                        treeInfo.right_children++;
                        right = new AABBNode(rightSide, rightSide_face_id, treeInfo, depth + 1);
                    }
                }
                else
                {
                    // should never happen
                    bMakeChildren = false;
                }
            }

            if (!bMakeChildren)
            {
                // Store the data directly if you want....
                for (int i = 0; i < vertexList.size(); i++)
                {
                    m_pVerts.emplace_back(vertexList[i]);
                }
                for (int i = 0; i < face_id_list.size(); i++)
                {
                    face_id.emplace_back(face_id_list[i]);
                }
            }
        }
    }
    ///////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    template <typename T>
    AABB_Tree<T>::AABB_Tree(const std::vector<Vector> &vertices_list)
    {
        treeInfo.max_tree_depth = 100;
        treeInfo.min_vertices = 72;
        treeInfo.curr_max_depth = 0;
        treeInfo.m_bPrune = true;
        std::vector<ptrdiff_t> face_id_list;
        face_id_list.reserve(vertices_list.size() / 3);
        for (int i = 0; i < vertices_list.size() / 3; i++)
            face_id_list.emplace_back(i);

        vertices_list_copy.assign(vertices_list.begin(), vertices_list.end());

        root = new aabb_impl::AABBNode<T>(vertices_list, face_id_list, treeInfo, 0);
        kdTree = 0;
        if (Enable_AABB_KDTree)
        {
            int seq;
            std::vector<Vector> v_list;
            std::vector<ptrdiff_t> local_face_id_list;
            v_list.clear();
            // sample some points
            Vector face_center;
            for (int i = 0; i < vertices_list.size() / 3; i++)
            {
                seq = rand() % 3;
                v_list.emplace_back(vertices_list[3 * i + seq]);
                local_face_id_list.emplace_back(i);
                face_center = (vertices_list[3 * i] + vertices_list[3 * i + 1] + vertices_list[3 * i + 2]) / (T)3;
                v_list.emplace_back(face_center);
                local_face_id_list.emplace_back(i);
            }
            kdTree = new KD_Tree<T, 3>(v_list, local_face_id_list);
        }
    }
    //////////////////////////////////////////////////////////////
    template <typename T>
    AABB_Tree<T>::~AABB_Tree()
    {
        if (root)
            delete root;
        if (kdTree)
            delete kdTree;
    }
    //////////////////////////////////////////////////////////////
    template <typename T>
    TinyVector<T, 3> AABB_Tree<T>::get_normal(ptrdiff_t face_ID)
    {
        return (vertices_list_copy[3 * face_ID + 1] - vertices_list_copy[3 * face_ID]).UnitCross(vertices_list_copy[3 * face_ID + 2] - vertices_list_copy[3 * face_ID]);
    }
    //////////////////////////////////////////////////////////////
    template <typename T>
    T AABB_Tree<T>::point_line_distance(const Vector &p, const Vector &a, const Vector &b, Vector &nearestP, T ref_dist)
    {
        // ref_dist is used for efficiency
        // reduce some computation
        Vector ab = b - a;
        Vector pa = a - p;
        Vector p_inline;
        T point_line_dist;
        T rst_dist;
        ab.Normalize();
        point_line_dist = (pa - pa.Dot(ab) * ab).Length();
        p_inline = p + pa - pa.Dot(ab) * ab;

        if (ref_dist >= 0 && point_line_dist > ref_dist)
        {
            // reduction of computation
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
    bool AABB_Tree<T>::inside_segment(const Vector &p, const Vector &a, const Vector &b)
    {
        Vector pa, pb;
        pa = a - p;
        pb = b - p;
        pa.Normalize();
        pb.Normalize();
        if (fabs(fabs(pa.Dot(pb)) - 1) > threshold)
            return false;
        if (pa.Dot(pb) < -0.5)
            return true;
        return false;
    }
    //////////////////////////////////////////////////////////////
    template <typename T>
    void AABB_Tree<T>::Traverse(aabb_impl::AABBNode<T> *node)
    {
        // statistic face number
        if (node == nullptr)
            return;
        if (node->left == nullptr && node->right == nullptr)
        {
            leafNodeNum++;
            assert(node->m_nNumVerts % 3 == 0);
            all_face_num += node->m_nNumVerts / 3;
            return;
        }
        if (node->left != nullptr)
            Traverse(node->left);
        if (node->right != nullptr)
            Traverse(node->right);
    }
    //////////////////////////////////////////////////////////////
    template <typename T>
    T AABB_Tree<T>::point_AABB_distance(const Vector &p, const aabb_impl::AABB<T> &box)
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

        Vector corner[8];
        if (x >= minx && x <= maxx && y >= miny && y <= maxy && z >= minz && z <= maxz)
        {
            // inside AABB
            return -1;
        }
        mindist = 0;
        if (x < minx)
        {
            mindist += (x - minx) * (x - minx);
        }
        else if (x > maxx)
        {
            mindist += (x - maxx) * (x - maxx);
        }

        if (y < miny)
        {
            mindist += (miny - y) * (miny - y);
        }
        else if (y > maxy)
        {
            mindist += (y - maxy) * (y - maxy);
        }

        if (z < minz)
        {
            mindist += (minz - z) * (minz - z);
        }
        else if (z > maxz)
        {
            mindist += (z - maxz) * (z - maxz);
        }

        return sqrt(mindist);
    }
    //////////////////////////////////////////////////////////////
    template <typename T>
    void AABB_Tree<T>::find_Nearest_Point(const Vector &p, Vector &nearestP, ptrdiff_t *face_ID)
    {
        T global_min_dist;
        Vector globalVector;
        ptrdiff_t globalFaceID;
        if (Enable_AABB_KDTree)
        {
            ptrdiff_t cur_face_id;
            kdTree->find_Nearest_Point(p, globalVector, &cur_face_id);
            globalFaceID = cur_face_id;
            global_min_dist = (globalVector - p).Length();
        }
        else
        {
            global_min_dist = (vertices_list_copy[0] - p).Length();
            globalVector = vertices_list_copy[0];
            globalFaceID = 0;
        }

        /*if ((last_rst - p).Length() < global_min_dist)
        {
            global_min_dist = (last_rst - p).Length();
            globalVector = last_rst;
            globalFaceID = last_rst_ID;
        }*/

        Traverse_Search(root, p, globalVector, global_min_dist, globalFaceID);
        nearestP = globalVector;
        if (face_ID != nullptr)
        {
            if (globalFaceID == -1)
            {
                // should never goes here
                assert(0);
                for (int i = 0; i < vertices_list_copy.size(); i++)
                {
                    if ((globalVector - vertices_list_copy[i]).Length() < 1e-10)
                    {
                        globalFaceID = i / 3;
                        break;
                    }
                }
            }
            *face_ID = globalFaceID;
        }
        // printf("minimum distance=%f\n", global_min_dist);
        // printf("nearest point=%f %f %f\n", nearestP(0), nearestP(1), nearestP(2));
        // last_rst = nearestP;
        // last_rst_ID = globalFaceID;
        // return global_min_dist;
    }
    //////////////////////////////////////////////////////////////
    template <typename T>
    void AABB_Tree<T>::find_Nearest_Point_by_traverse(const Vector &p, Vector &nearestP, ptrdiff_t *face_ID)
    {
        // this function calculate distance by calculate distance from query to every triangle faces
        // just for test
        T global_min_dist;
        Vector globalVector;
        int globalFaceID;
        global_min_dist = (vertices_list_copy[0] - p).Length();
        globalVector = vertices_list_copy[0];
        globalFaceID = 0;

        Vector local_nearest;
        T local_dist;
        for (int i = 0; i < vertices_list_copy.size() / 3; i++)
        {
            local_dist = point_tri_distance_refine(p, vertices_list_copy[3 * i], vertices_list_copy[3 * i + 1], vertices_list_copy[3 * i + 2], local_nearest, -1);
            if (local_dist < global_min_dist)
            {
                global_min_dist = local_dist;
                globalVector = local_nearest;
                globalFaceID = i;
            }
        }

        nearestP = globalVector;
        if (face_ID != nullptr)
        {
            if (globalFaceID == -1)
            {
                assert(0);
                // should never goes here
                for (int i = 0; i < vertices_list_copy.size(); i++)
                {
                    if ((globalVector - vertices_list_copy[i]).Length() < 1e-10)
                    {
                        globalFaceID = i / 3;
                        break;
                    }
                }
            }
            *face_ID = globalFaceID;
        }
        // return global_min_dist;
    }
    //////////////////////////////////////////////////////////////
    template <typename T>
    void AABB_Tree<T>::Traverse_Search(aabb_impl::AABBNode<T> *node, const Vector &p, Vector &globalVector, T &global_min_dist, ptrdiff_t &globalFaceID)
    {
        if (node == nullptr)
            return;
        if (point_AABB_distance(p, node->m_box) > global_min_dist)
            return;
        if (node->left == nullptr && node->right == nullptr)
        {
            // search entire node
            search_entire_node(node, p, globalVector, global_min_dist, globalFaceID);
        }
        if (node->left != nullptr)
            Traverse_Search(node->left, p, globalVector, global_min_dist, globalFaceID);
        if (node->right != nullptr)
            Traverse_Search(node->right, p, globalVector, global_min_dist, globalFaceID);
    }
    //////////////////////////////////////////////////////////////
    template <typename T>
    T AABB_Tree<T>::search_entire_node(aabb_impl::AABBNode<T> *node, const Vector &p, Vector &globalVector, T &global_min_dist, ptrdiff_t &globalFaceID)
    {
        T min_dist, temp_rst;
        ptrdiff_t cur_face_id;
        Vector min_point;
        Vector rst;
        if (node == nullptr)
            return -1;
        if (node->m_nNumVerts < 3)
            return -1;

        min_dist = point_tri_distance_refine(p, node->m_pVerts[0], node->m_pVerts[1], node->m_pVerts[2], min_point, -1);
        cur_face_id = node->face_id[0];

        for (int i = 1; i < node->m_nNumVerts / 3; i++)
        {
            temp_rst = point_tri_distance_refine(p, node->m_pVerts[3 * i], node->m_pVerts[3 * i + 1], node->m_pVerts[3 * i + 2], rst, min_dist);
            if (temp_rst < min_dist)
            {
                min_dist = temp_rst;
                min_point = rst;
                cur_face_id = node->face_id[i];
            }
        }
        if (min_dist < global_min_dist)
        {
            global_min_dist = min_dist;
            globalVector = min_point;
            globalFaceID = cur_face_id;
        }
        return global_min_dist;
    }
    //////////////////////////////////////////////////////////////
    template <typename T>
    T AABB_Tree<T>::point_tri_distance_refine(const Vector &p, const Vector &a, const Vector &b, const Vector &c, Vector &nearestP, T ref_dist)
    {
        Vector center = (a + b + c) / (T)3;
        Vector ab, ac, bc;
        Vector nab, nbc, nac;
        Vector tri_face_normal, temp, p_projected;
        T nearestDist;
        T dist_to_plane;
        int state_code;
        ab = b - a;
        ac = c - a;
        bc = c - b;
        ab.Normalize();
        ac.Normalize();
        bc.Normalize();
        tri_face_normal = ab.Cross(ac);
        tri_face_normal.Normalize();

        temp = a - p;
        dist_to_plane = fabs(temp.Dot(tri_face_normal));
        if (ref_dist >= 0 && ref_dist < dist_to_plane)
        {
            // reduction for efficiency
            return dist_to_plane;
        }

        nab = center - a;
        nbc = center - b;
        nac = center - a;
        nab = nab - nab.Dot(ab) * ab;
        nbc = nbc - nbc.Dot(bc) * bc;
        nac = nac - nac.Dot(ac) * ac;
        nab.Normalize();
        nbc.Normalize();
        nac.Normalize();
        p_projected = p + temp.Dot(tri_face_normal) * tri_face_normal;

        T d1, d2, d3, d4, d5;
        d1 = tri_face_normal.Dot(ab);
        d2 = tri_face_normal.Dot(ac);
        d3 = tri_face_normal.Dot(bc);
        d4 = tri_face_normal.Dot(p_projected - c);
        d5 = tri_face_normal.Dot(p_projected - b);

        state_code = 0;
        if ((p_projected - a).Dot(nab) > threshold)
            state_code += 1;
        if ((p_projected - b).Dot(nbc) > threshold)
            state_code += 2;
        if ((p_projected - c).Dot(nac) > threshold)
            state_code += 4;
        if (state_code == 7)
        {
            // inside triangle
            nearestP = p_projected;
            nearestDist = dist_to_plane;
        }
        else
        {
            // out of triangle
            // find nearest point on edges or vertices
            Vector tempP;
            /*switch (state_code)
            {
            case 1:
                nearestDist = point_line_distance(p, a, b, tempP, -1);
                nearestP = tempP;
                break;
            case 2:
                nearestDist = point_line_distance(p, b, c, tempP, -1);
                nearestP = tempP;
                break;
            case 3:
                nearestDist = (p - b).Length();
                nearestP = b;
                break;
            case 4:
                nearestDist = point_line_distance(p, a, c, tempP, -1);
                nearestP = tempP;
                break;
            case 5:
                nearestDist = (p - a).Length();
                nearestP = a;
                break;
            case 6:
                nearestDist = (p - c).Length();
                nearestP = c;
                break;
            default:
                //never goes here
                break;
            }*/
            T local_dist;
            Vector local_rst;
            switch (state_code)
            {
            case 1:
                // check ac bc
                nearestDist = point_line_distance(p, a, c, tempP, -1);
                nearestP = tempP;

                local_dist = point_line_distance(p, b, c, local_rst, nearestDist);
                if (local_dist < nearestDist)
                {
                    nearestDist = local_dist;
                    nearestP = local_rst;
                }
                break;
            case 2:
                // check ab ac
                nearestDist = point_line_distance(p, a, c, tempP, -1);
                nearestP = tempP;

                local_dist = point_line_distance(p, a, b, local_rst, nearestDist);
                if (local_dist < nearestDist)
                {
                    nearestDist = local_dist;
                    nearestP = local_rst;
                }
                break;
            case 3:
                // check ac
                nearestDist = point_line_distance(p, a, c, tempP, -1);
                nearestP = tempP;
                break;
            case 4:
                // check bc ab
                nearestDist = point_line_distance(p, b, c, tempP, -1);
                nearestP = tempP;

                local_dist = point_line_distance(p, a, b, local_rst, nearestDist);
                if (local_dist < nearestDist)
                {
                    nearestDist = local_dist;
                    nearestP = local_rst;
                }
                break;
            case 5:
                // check bc
                nearestDist = point_line_distance(p, b, c, tempP, -1);
                nearestP = tempP;
                break;
            case 6:
                // check ab
                nearestDist = point_line_distance(p, a, b, tempP, -1);
                nearestP = tempP;
                break;
            default:
                // never goes here
                // check ab ac bc
                nearestDist = point_line_distance(p, b, c, tempP, -1);
                nearestP = tempP;

                local_dist = point_line_distance(p, a, b, local_rst, nearestDist);
                if (local_dist < nearestDist)
                {
                    nearestDist = local_dist;
                    nearestP = local_rst;
                }

                local_dist = point_line_distance(p, a, c, local_rst, nearestDist);
                if (local_dist < nearestDist)
                {
                    nearestDist = local_dist;
                    nearestP = local_rst;
                }
                break;
            }
        }
        return nearestDist;
    }
    //////////////////////////////////////////////////////////////
}
