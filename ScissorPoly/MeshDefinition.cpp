#include "MeshDefinition.h"
#include "Eigen/Dense"

OpenVolumeMesh::HalfFaceHandle oppositeHalffaceInCell(OpenVolumeMesh::HalfFaceHandle hf, VolumeMesh* mesh_)
{
	OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face(mesh_->face_handle(hf));
	std::vector<OpenVolumeMesh::HalfEdgeHandle> hfe_Vec = face.halfedges();
	std::vector<int> edge_id;
	for(unsigned int i=0;i<hfe_Vec.size();++i)
	{
		edge_id.push_back( mesh_->edge_handle(hfe_Vec[i]).idx() );
	}

	OpenVolumeMesh::CellHandle ch = mesh_->incident_cell(hf);
	OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hff_vec = cell.halffaces();
	OpenVolumeMesh::HalfFaceHandle result;
	for(unsigned int i=0;i<hff_vec.size();++i)
	{
		face = mesh_->face(mesh_->face_handle(hff_vec[i]));
		hfe_Vec = face.halfedges();
		int count = 0;
		for(unsigned int j=0;j<hfe_Vec.size();++j)
		{
			if( std::find( edge_id.begin(),edge_id.end(), mesh_->edge_handle(hfe_Vec[j]).idx()) != edge_id.end()  )
			{
				++count;
			}
			if(count > 0)
			{
				break;
			}
		}
		if(0 == count)
		{
			result = hff_vec[i]; break;
		}
	}
	return result;
}

OpenVolumeMesh::EdgeHandle oppositeEdgeInFace(OpenVolumeMesh::HalfFaceHandle hf, OpenVolumeMesh::EdgeHandle e, VolumeMesh* mesh_)
{
	OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face( mesh_->face_handle(hf) );
	std::vector<OpenVolumeMesh::HalfEdgeHandle> heh_Vec = face.halfedges();
	OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge(e);
	for(unsigned int i=0;i<heh_Vec.size();++i)
	{
		OpenVolumeMesh::OpenVolumeMeshEdge edge2 = mesh_->edge(mesh_->edge_handle(heh_Vec[i]));
		if(edge2.from_vertex() != edge.from_vertex() && edge2.from_vertex() != edge.to_vertex() 
			&& edge2.to_vertex() != edge.from_vertex() && edge2.to_vertex() != edge.to_vertex() )
		{
			return mesh_->edge_handle(heh_Vec[i]);
		}
	}
	return NULL;
}

OpenVolumeMesh::EdgeHandle oppositeEdgeInOppositeFace(OpenVolumeMesh::HalfFaceHandle hf,OpenVolumeMesh::EdgeHandle e, VolumeMesh* mesh_)
{
	OpenVolumeMesh::CellHandle ch = mesh_->incident_cell(hf);
	OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hff_vec = cell.halffaces();
	OpenVolumeMesh::HalfFaceHandle temp_hf;
	for( unsigned int i=0; i<hff_vec.size(); ++i )
	{
		if(hff_vec[i] == hf) continue;

		OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face( mesh_->face_handle(hff_vec[i]) );
		std::vector<OpenVolumeMesh::HalfEdgeHandle> heh_Vec = face.halfedges();
		bool edge_OK = false;
		for(unsigned int j=0;j<heh_Vec.size();++j)
		{
			if(mesh_->edge_handle(heh_Vec[j]) == e)
			{
				temp_hf = hff_vec[i];
				edge_OK = true;
				break;
			}
		}
		if(edge_OK)
		{break;}
	}

	return oppositeEdgeInFace(temp_hf,e,mesh_);
}



std::vector<OpenVolumeMesh::HalfFaceHandle> adjHalfFaceinSheetinCell( OpenVolumeMesh::HalfFaceHandle hf, OpenVolumeMesh::EdgeHandle eh[2] , VolumeMesh* mesh_)
{
	int edgeId[2] = {eh[0].idx(), eh[1].idx()};
	std::vector<OpenVolumeMesh::HalfFaceHandle> result;
	OpenVolumeMesh::CellHandle ch = mesh_->incident_cell(hf);
	if(!ch.is_valid())
	{
		return result;
	}
	result.resize(3);
	OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hff_vec = cell.halffaces();


	OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face(mesh_->face_handle(hf));
	std::vector<OpenVolumeMesh::HalfEdgeHandle> hfe_Vec = face.halfedges();
	std::vector<int> hf_edgeID;
	for(unsigned int j=0;j<hfe_Vec.size();++j)
	{
		hf_edgeID.push_back( mesh_->edge_handle(hfe_Vec[j]).idx() );
	}
	int temp_edgeId = 0;
	int count = 0;
	int adjHF = 0;
	for(unsigned int i=0;i<hff_vec.size();++i)
	{
		face = mesh_->face(mesh_->face_handle(hff_vec[i]));
		hfe_Vec = face.halfedges();
		count = 0 ;
		for(unsigned int j=0;j<hfe_Vec.size();++j)
		{
			temp_edgeId = mesh_->edge_handle(hfe_Vec[j]).idx();
			if(temp_edgeId == edgeId[0] || temp_edgeId == edgeId[1])
			{
				++count;
			}
		}
		if(count==1)
		{
			result[adjHF] = hff_vec[i]; ++adjHF;
		}
		else if(count == 0)
		{
			//check the opposite edge
			for(unsigned int j=0;j<hfe_Vec.size();++j)
			{
				temp_edgeId = mesh_->edge_handle(hfe_Vec[j]).idx();
				for(unsigned int k=0;k<hf_edgeID.size();++k)
				{
					if(temp_edgeId == hf_edgeID[k])
					{
						++count;
						break;
					}
				}
				if(count >0 )
				{
					break;
				}
			}
			if(count == 0)
			{
				result[2] = hff_vec[i];
			}
		}
	}

	return result;
}

std::vector<OpenVolumeMesh::HalfFaceHandle > threeOppositeFacesinCell(OpenVolumeMesh::HalfFaceHandle hf, VolumeMesh* mesh_)
{
	OpenVolumeMesh::CellHandle ch = mesh_->incident_cell(hf);
	OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hff_Vec = cell.halffaces();
	std::vector<OpenVolumeMesh::HalfFaceHandle> result;
	std::vector<int> temp_result;
	result.push_back(hff_Vec[0]);
	temp_result.push_back(result[0].idx()); temp_result.push_back(oppositeHalffaceInCell(result[0],mesh_).idx());
	int temp_id = 0; bool exist = false;
	for(unsigned int i=1;i<hff_Vec.size();++i)
	{
		temp_id = hff_Vec[i].idx();
		exist = false;
		for(unsigned int j=0;j<temp_result.size();++j)
		{
			if(temp_id == temp_result[j])
			{
				exist = true;
				break;
			}
		}
		if( !exist )
		{
			result.push_back(hff_Vec[i]);
			temp_result.push_back(hff_Vec[i].idx()); temp_result.push_back(oppositeHalffaceInCell(hff_Vec[i],mesh_).idx());
		}
	}

	return result;
}


OpenVolumeMesh::EdgeHandle edgeBetweenTwoHalfFace(OpenVolumeMesh::HalfFaceHandle hf[2], VolumeMesh* mesh_)
{
	OpenVolumeMesh::OpenVolumeMeshFace face0 = mesh_->face(mesh_->face_handle(hf[0]));
	std::vector<OpenVolumeMesh::HalfEdgeHandle> hfe_Vec0 = face0.halfedges();

	OpenVolumeMesh::OpenVolumeMeshFace face1 = mesh_->face(mesh_->face_handle(hf[1]));
	std::vector<OpenVolumeMesh::HalfEdgeHandle> hfe_Vec1 = face1.halfedges();

	OpenVolumeMesh::EdgeHandle result;
	for(unsigned int i=0;i<hfe_Vec0.size();++i)
	{
		for(unsigned int j=0;j<hfe_Vec1.size();++j)
		{
			if( mesh_->edge_handle(hfe_Vec0[i]).idx()  == mesh_->edge_handle(hfe_Vec1[j]).idx() )
			{
				result = mesh_->edge_handle(hfe_Vec0[i]);
				break;
			}
		}
	}
	return result;
}

OpenVolumeMesh::EdgeHandle edgeBetweenTwoHalfFace(OpenVolumeMesh::HalfFaceHandle hf1,OpenVolumeMesh::HalfFaceHandle hf2, VolumeMesh* mesh_)
{
	OpenVolumeMesh::OpenVolumeMeshFace face0 = mesh_->face(mesh_->face_handle(hf1));
	std::vector<OpenVolumeMesh::HalfEdgeHandle> hfe_Vec0 = face0.halfedges();

	OpenVolumeMesh::OpenVolumeMeshFace face1 = mesh_->face(mesh_->face_handle(hf2));
	std::vector<OpenVolumeMesh::HalfEdgeHandle> hfe_Vec1 = face1.halfedges();

	OpenVolumeMesh::EdgeHandle result;
	for(unsigned int i=0;i<hfe_Vec0.size();++i)
	{
		for(unsigned int j=0;j<hfe_Vec1.size();++j)
		{
			if( mesh_->edge_handle(hfe_Vec0[i]).idx()  == mesh_->edge_handle(hfe_Vec1[j]).idx() )
			{
				result = mesh_->edge_handle(hfe_Vec0[i]);
				break;
			}
		}
	}
	return result;
}

std::vector<std::vector<OpenVolumeMesh::HalfFaceHandle > > HalfFaceinSheet(OpenVolumeMesh::HalfFaceHandle hf, VolumeMesh* mesh_)
{
	OpenVolumeMesh::CellHandle ch = mesh_->incident_cell(hf);
	OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hff_vec = cell.halffaces();

	OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face(mesh_->face_handle(hf));
	std::vector<OpenVolumeMesh::HalfEdgeHandle> hfe_Vec = face.halfedges();

	std::vector<int> edge_id;
	for(unsigned int i=0;i<hfe_Vec.size();++i)
	{
		edge_id.push_back( mesh_->edge_handle(hfe_Vec[i]).idx() );
	}
	int count =0; int temp_edge_id=0;
	OpenVolumeMesh::HalfFaceHandle oppositeHF;
	std::vector< OpenVolumeMesh::HalfFaceHandle > adjHF;
	for(unsigned int i=0; i<hff_vec.size(); ++i )
	{
		face = mesh_->face(mesh_->face_handle(hff_vec[i]));
		hfe_Vec = face.halfedges();
		count = 0;
		for(unsigned int j=0; j<hfe_Vec.size(); ++j)
		{
			temp_edge_id = mesh_->edge_handle(hfe_Vec[j]).idx();
			for(unsigned int k=0;k<edge_id.size();++k)
			{
				if(edge_id[k] == temp_edge_id)
				{
					++count; break;
				}
			}
			if(count > 1)
			{
				break;
			}
		}

		if(count == 0)//opposite
		{
			oppositeHF = hff_vec[i];
		}
		else if(count == 1)//adjacent
		{
			adjHF.push_back(hff_vec[i]);
		}
		else if(count >1 )//self
		{

		}
	}

	adjustOrderofFourAdjHfinCell(adjHF,mesh_);

	std::vector<std::vector<OpenVolumeMesh::HalfFaceHandle > > result;
	result.push_back(adjHF);

	std::vector<OpenVolumeMesh::HalfFaceHandle > other; other.resize(4);
	other[0] = hf; other[2] = oppositeHF; 
	other[1] = adjHF[0]; other[3] = adjHF[2]; 
	result.push_back(other);
	other[1] = adjHF[1]; other[3] = adjHF[3];
	result.push_back(other);

	return result;
}

std::vector<OpenVolumeMesh::HalfFaceHandle > FourHalffacesinCell(OpenVolumeMesh::CellHandle ch, std::vector<OpenVolumeMesh::EdgeHandle> eh_Vec, VolumeMesh* mesh_)
{
	OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hff_Vec = cell.halffaces();
	OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face( mesh_->face_handle(hff_Vec[0]) );
	std::vector<OpenVolumeMesh::HalfEdgeHandle> heh_Vec = face.halfedges();
	std::vector<OpenVolumeMesh::HalfFaceHandle> result;
	for(unsigned int i=0;i<hff_Vec.size();++i)
	{
		face = mesh_->face( mesh_->face_handle(hff_Vec[i]) );
		heh_Vec = face.halfedges();
		bool added = false;
		for( unsigned int j=0; j<heh_Vec.size(); ++j )
		{
			int temp_eh_ID = mesh_->edge_handle(heh_Vec[j]).idx();
			for( unsigned int k=0; k<eh_Vec.size(); ++k)
			{
				if( eh_Vec[k] == temp_eh_ID )
				{
					result.push_back(hff_Vec[i]);
					added = true;
					break;
				}
			}
			if( added )
			{
				break;
			}
		}
	}

	adjustOrderofFourAdjHfinCell(result,mesh_);

	return result;
}

void adjustOrderofFourAdjHfinCell(std::vector<OpenVolumeMesh::HalfFaceHandle>& adjHF, VolumeMesh* mesh_)
{
	OpenVolumeMesh::HalfFaceHandle oppositeHF = oppositeHalffaceInCell(adjHF[0],mesh_);
	for(unsigned i=1;i<adjHF.size();++i)
	{
		if( oppositeHF.idx() == adjHF[i].idx() )
		{
			if( i != 2)
			{
				adjHF[i] = adjHF[2];
				adjHF[2] = oppositeHF;
			}
			break;
		}
	}
}

std::vector<OpenVolumeMesh::HalfFaceHandle > AdjFourHalffacesinCell(OpenVolumeMesh::HalfFaceHandle hf, VolumeMesh* mesh_)
{
	OpenVolumeMesh::CellHandle ch = mesh_->incident_cell(hf);
	OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hff_Vec = cell.halffaces();
	OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face( mesh_->face_handle(hf) );
	std::vector<OpenVolumeMesh::HalfEdgeHandle> heh_Vec = face.halfedges();
	std::vector<int> edge_id(4);
	for(unsigned int i=0;i<heh_Vec.size();++i)
	{
		edge_id[i] = mesh_->edge_handle(heh_Vec[i]).idx();
	}

	std::vector<OpenVolumeMesh::HalfFaceHandle> result;
	std::vector<bool> HFVisited(hff_Vec.size(),false);
	for(unsigned int j=0;j<edge_id.size();++j)
	{
		for(unsigned int i=0;i<hff_Vec.size();++i)
		{
			if( hff_Vec[i].idx() != hf.idx() && !HFVisited[i] )
			{
				face = mesh_->face( mesh_->face_handle(hff_Vec[i]) );
				heh_Vec = face.halfedges();
				for( unsigned int k=0; k<heh_Vec.size(); ++k )
				{
					if( mesh_->edge_handle(heh_Vec[k]).idx() == edge_id[j] )
					{
						result.push_back(hff_Vec[i]);
						HFVisited[i] = true;
						break;
					}
				}
			}
		}
	}

	return result;
}


std::vector<OpenVolumeMesh::Geometry::Vec3d> vertexOnFace(OpenVolumeMesh::HalfFaceHandle hf, VolumeMesh* mesh_)
{
	OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hf);
	std::vector<OpenVolumeMesh::Geometry::Vec3d> result;
	for( hfv_it; hfv_it; ++hfv_it )
	{
		result.push_back( mesh_->vertex(*hfv_it) );
	}
	return result;
}

std::vector<OpenVolumeMesh::VertexHandle> VertexOnFaceinOrder(OpenVolumeMesh::HalfFaceHandle hf, VolumeMesh* mesh_)
{
	std::vector<OpenVolumeMesh::VertexHandle> result;
	return result;
}

std::vector<int> edgeOnFace(OpenVolumeMesh::HalfFaceHandle hf, VolumeMesh* mesh_)
{
	OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face(mesh_->face_handle(hf));
	std::vector<OpenVolumeMesh::HalfEdgeHandle> hfe_Vec = face.halfedges();

	std::vector<int> result(hfe_Vec.size());
	for(unsigned int i=0;i<hfe_Vec.size();++i)
	{
		result[i] = mesh_->edge_handle( hfe_Vec[i]).idx();
	}
	return result;
}

std::vector<OpenVolumeMesh::VertexHandle> vertexOnEdge(OpenVolumeMesh::EdgeHandle eh, VolumeMesh* mesh_)
{
	std::vector<OpenVolumeMesh::VertexHandle> result;
	OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge(eh);
	result.push_back( edge.from_vertex() );
	result.push_back( edge.to_vertex() );
	return result;
}


OpenVolumeMesh::Geometry::Vec3d midPointOfEdge(OpenVolumeMesh::EdgeHandle eh, VolumeMesh* mesh_)
{
	OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge(eh);
	OpenVolumeMesh::VertexHandle vh1 = edge.from_vertex();
	OpenVolumeMesh::VertexHandle vh2 = edge.to_vertex();

	return (mesh_->vertex(vh1) + mesh_->vertex(vh2))*0.5;
}

OpenVolumeMesh::VertexHandle adjVertexBetweenTwoFace(OpenVolumeMesh::HalfFaceHandle hf1,OpenVolumeMesh::VertexHandle v1,OpenVolumeMesh::HalfFaceHandle hf2, VolumeMesh* mesh_)
{
	assert( mesh_->incident_cell(hf1) == mesh_->incident_cell(hf2) );

	OpenVolumeMesh::CellHandle ch = mesh_->incident_cell(hf1);
	OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hff_Vec = cell.halffaces();
	OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hff_Vec[0]);
	OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face( mesh_->face_handle(hff_Vec[0]) );
	std::vector<OpenVolumeMesh::HalfEdgeHandle> hfe_Vec = face.halfedges();
	std::vector<OpenVolumeMesh::VertexHandle> twoVertexOnEdge;
	bool vertexOnFace = false;
	for(unsigned int i=0; i<hff_Vec.size(); ++i)
	{
		if(hff_Vec[i].idx() != hf1.idx() && hff_Vec[i].idx() != hf2.idx())
		{ 
			hfv_it = mesh_->hfv_iter(hff_Vec[i]);
			for( hfv_it; hfv_it; ++hfv_it )
			{
				if(hfv_it->idx() == v1.idx())
				{
					vertexOnFace = true;
					break;
				}
			}
			if(vertexOnFace)
			{
				face = mesh_->face( mesh_->face_handle(hff_Vec[i]) );
				hfe_Vec = face.halfedges();
				for(unsigned int j=0;j<hfe_Vec.size();++j)
				{
					twoVertexOnEdge = vertexOnEdge(mesh_->edge_handle(hfe_Vec[j]),mesh_);
					if(twoVertexOnEdge[0].idx() == v1.idx())
					{
						hfv_it = mesh_->hfv_iter(hf2);
						for( hfv_it; hfv_it; ++hfv_it )
						{
							if(hfv_it->idx() == twoVertexOnEdge[1].idx())
							{
								return twoVertexOnEdge[1];
							}
						}

					}
					else if(twoVertexOnEdge[1].idx() == v1.idx())
					{
						hfv_it = mesh_->hfv_iter(hf2);
						for( hfv_it; hfv_it; ++hfv_it )
						{
							if(hfv_it->idx() == twoVertexOnEdge[0].idx())
							{
								return twoVertexOnEdge[0];
							}
						}
					}
				}
			}
		}
	}

	return twoVertexOnEdge[0];//failed
}

std::vector<OpenVolumeMesh::HalfFaceHandle> adjTwoHalfFaceInCell(std::vector<OpenVolumeMesh::HalfFaceHandle> hf2,std::vector<OpenVolumeMesh::VertexHandle> v2, VolumeMesh* mesh_)
{
	assert( hf2.size() == 2 && v2.size() == 2 );
	OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hf2[0]);
	std::vector<OpenVolumeMesh::HalfFaceHandle> result(2);
	OpenVolumeMesh::CellHandle ch = mesh_->incident_cell(hf2[0]);
	OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(mesh_->incident_cell(hf2[0]));
	std::vector<OpenVolumeMesh::HalfFaceHandle> hff_Vec = cell.halffaces();
	int vertexOnFaceCount =0;
	for(unsigned int i=0;i<hff_Vec.size();++i)
	{
		if(hff_Vec[i].idx() != hf2[0].idx() && hff_Vec[i].idx() != hf2[1].idx() )
		{
			hfv_it = mesh_->hfv_iter(hff_Vec[i]);
			vertexOnFaceCount = 0;
			for(hfv_it;hfv_it;++hfv_it)
			{
				if( hfv_it->idx() == v2[0].idx() )//keep  the order between vertex and hf
				{
					result[0] =hff_Vec[i];
					break;
				}
				else if( hfv_it->idx() == v2[1].idx())
				{
					result[1] =hff_Vec[i];
					break;
				}
			}
		}
	}
	return result;
}

VertexEdgeOnDoublet vertexOnOppositeFace(OpenVolumeMesh::HalfFaceHandle hf, OpenVolumeMesh::VertexHandle v, VolumeMesh* mesh_)
{
	OpenVolumeMesh::CellHandle ch =  mesh_->incident_cell(hf);
	OpenVolumeMesh::OpenVolumeMeshCell cell =mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hff_Vec = cell.halffaces();
	OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hf);
	OpenVolumeMesh::HalfFaceHandle temp_hf; 
	bool haveFace = false;
	for(unsigned int i=0;i<hff_Vec.size();++i)
	{
		if(hff_Vec[i].idx() == hf.idx() )
		{
			continue;
		}

		hfv_it = mesh_->hfv_iter(hff_Vec[i]); 
		for( hfv_it; hfv_it; ++hfv_it )
		{
			if(hfv_it->idx() == v.idx() )
			{
				haveFace = true;
				break;
			}
		}
		if( haveFace )
		{
			temp_hf = hff_Vec[i];
			break;
		}
	}

	OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face(mesh_->face_handle(temp_hf));
	std::vector<OpenVolumeMesh::HalfEdgeHandle> hfe_Vec = face.halfedges();
	OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge(mesh_->edge_handle(hfe_Vec[0]));
	std::vector<OpenVolumeMesh::EdgeHandle> temp_eh_Vec;
	std::vector<OpenVolumeMesh::VertexHandle> temp_vh_Vec;
	OpenVolumeMesh::EdgeHandle temp_eh;
	for( unsigned int i=0; i<hfe_Vec.size(); ++i )
	{
		temp_eh = mesh_->edge_handle(hfe_Vec[i]);
		edge = mesh_->edge(temp_eh);
		if(edge.from_vertex() == v )
		{
			temp_eh_Vec.push_back( temp_eh );
			temp_vh_Vec.push_back(edge.to_vertex());
		}
		else if(edge.to_vertex() == v)
		{
			temp_eh_Vec.push_back( temp_eh );
			temp_vh_Vec.push_back(edge.from_vertex());
		}
	}

	VertexEdgeOnDoublet result;
	hfv_it = mesh_->hfv_iter(hf); 
	for( hfv_it; hfv_it; ++hfv_it )
	{
		if(hfv_it->idx() == temp_vh_Vec[0].idx())
		{
			result = VertexEdgeOnDoublet(temp_vh_Vec[1],temp_eh_Vec[1]);
			break;
		}
		else if( hfv_it->idx() == temp_vh_Vec[1].idx() )
		{
			result = VertexEdgeOnDoublet(temp_vh_Vec[0],temp_eh_Vec[0]);
			break;
		}
	}
	return result;

}

OpenVolumeMesh::VertexHandle VertexHandleOnOppositeFace(OpenVolumeMesh::HalfFaceHandle hf, OpenVolumeMesh::VertexHandle v, VolumeMesh* mesh_)
{
	OpenVolumeMesh::CellHandle ch =  mesh_->incident_cell(hf);
	OpenVolumeMesh::OpenVolumeMeshCell cell =mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hff_Vec = cell.halffaces();
	OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hf);
	OpenVolumeMesh::HalfFaceHandle temp_hf; 
	bool haveFace = false;
	for(unsigned int i=0;i<hff_Vec.size();++i)
	{
		if(hff_Vec[i].idx() == hf.idx() )
		{
			continue;
		}

		hfv_it = mesh_->hfv_iter(hff_Vec[i]); 
		for( hfv_it; hfv_it; ++hfv_it )
		{
			if(hfv_it->idx() == v.idx() )
			{
				haveFace = true;
				break;
			}
		}
		if( haveFace )
		{
			temp_hf = hff_Vec[i];
			break;
		}
	}

	OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face(mesh_->face_handle(temp_hf));
	std::vector<OpenVolumeMesh::HalfEdgeHandle> hfe_Vec = face.halfedges();
	OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge(mesh_->edge_handle(hfe_Vec[0]));
	std::vector<OpenVolumeMesh::EdgeHandle> temp_eh_Vec;
	std::vector<OpenVolumeMesh::VertexHandle> temp_vh_Vec;
	OpenVolumeMesh::EdgeHandle temp_eh;
	for( unsigned int i=0; i<hfe_Vec.size(); ++i )
	{
		temp_eh = mesh_->edge_handle(hfe_Vec[i]);
		edge = mesh_->edge(temp_eh);
		if(edge.from_vertex() == v )
		{
			temp_eh_Vec.push_back( temp_eh );
			temp_vh_Vec.push_back(edge.to_vertex());
		}
		else if(edge.to_vertex() == v)
		{
			temp_eh_Vec.push_back( temp_eh );
			temp_vh_Vec.push_back(edge.from_vertex());
		}
	}

	OpenVolumeMesh::VertexHandle result;
	hfv_it = mesh_->hfv_iter(hf); 
	for( hfv_it; hfv_it; ++hfv_it )
	{
		if(hfv_it->idx() == temp_vh_Vec[0].idx())
		{
			result = temp_vh_Vec[1];
			break;
		}
		else if( hfv_it->idx() == temp_vh_Vec[1].idx() )
		{
			result = temp_vh_Vec[0];
			break;
		}
	}

	return result;
}

OpenVolumeMesh::VertexHandle oppositeVertexOnFace(OpenVolumeMesh::HalfFaceHandle hf,OpenVolumeMesh::VertexHandle v, VolumeMesh* mesh_)
{
	OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hf);
	OpenVolumeMesh::VertexHandle result;
	std::vector<OpenVolumeMesh::VertexHandle> v_Vec;
	int index =0; int count =0;
	for(hfv_it;hfv_it;++hfv_it)
	{
		v_Vec.push_back(*hfv_it);
		if( hfv_it->idx() == v.idx() )
		{
			index = count;
		}
		++count;
	}
	result = v_Vec[(index +2)%4];
	return result;
}

OpenVolumeMesh::HalfFaceHandle findHFWithFourVertex(OpenVolumeMesh::HalfFaceHandle hf, std::vector<OpenVolumeMesh::VertexHandle> v_Vec, VolumeMesh* mesh_)
{
	OpenVolumeMesh::CellHandle ch =  mesh_->incident_cell(hf);
	OpenVolumeMesh::OpenVolumeMeshCell cell =mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hff_Vec = cell.halffaces();
	OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hf);
	OpenVolumeMesh::HalfFaceHandle result;
	for(unsigned int i=0;i<hff_Vec.size();++i)
	{
		hfv_it = mesh_->hfv_iter(hff_Vec[i]);
		int samNum = 0;
		for( hfv_it; hfv_it; ++hfv_it )
		{
			for(unsigned int j=0;j<v_Vec.size();++j)
			{
				if(hfv_it->idx() == v_Vec[j])
				{
					++samNum;
					break;
				}
			}
			if(samNum == v_Vec.size())
			{
				return hff_Vec[i];
			}
		}

	}
	return result;
}

OpenVolumeMesh::HalfFaceHandle HFonAdjCell(OpenVolumeMesh::HalfFaceHandle hf, std::vector<OpenVolumeMesh::HalfFaceHandle> TOHF, VolumeMesh* mesh_)
{
	assert(TOHF.size() == 2);
	assert(mesh_->incident_cell(TOHF[0]).idx() == mesh_->incident_cell(TOHF[1]).idx());

	OpenVolumeMesh::CellHandle ch = mesh_->incident_cell(TOHF[0]);
	OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hff_Vec = cell.halffaces();
	OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hf);
	std::vector<int> v_Vec;
	for(hfv_it;hfv_it;++hfv_it)
	{
		v_Vec.push_back(hfv_it->idx());
	}
	int vCount =0;
	OpenVolumeMesh::HalfFaceHandle result;
	for(unsigned int i=0;i<hff_Vec.size();++i)
	{
		if(hff_Vec[i].idx() != TOHF[0].idx() && hff_Vec[i].idx() != TOHF[1].idx())
		{
			hfv_it = mesh_->hfv_iter(hff_Vec[i]);
			vCount = 0;
			for(hfv_it;hfv_it;++hfv_it)
			{
				for(unsigned int j=0;j<v_Vec.size();++j)
				{
					if(hfv_it->idx() == v_Vec[j])
					{
						++vCount;
						break;
					}
				}
			}
			if(vCount == 2)
			{
				result = hff_Vec[i];
				break;
			}
		}
	}
	return result;
}

std::vector<twoOppositeVertex> twoPairOppositeVertexOnTwoFace(OpenVolumeMesh::HalfFaceHandle hf1,OpenVolumeMesh::HalfFaceHandle hf2,OpenVolumeMesh::EdgeHandle eh,OpenVolumeMesh::VertexHandle v1,OpenVolumeMesh::VertexHandle v2, VolumeMesh* mesh_)
{
	std::vector<twoOppositeVertex> result(2);

	OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face(mesh_->face_handle(hf1));
	std::vector<OpenVolumeMesh::HalfEdgeHandle> heh_Vec = face.halfedges();
	OpenVolumeMesh::EdgeHandle temp_eh = mesh_->edge_handle(heh_Vec[0]);
	OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge(temp_eh);
	for( unsigned int i=0; i<heh_Vec.size(); ++i )
	{
		temp_eh = mesh_->edge_handle(heh_Vec[i]);
		if( temp_eh.idx() != eh.idx() )
		{
			edge = mesh_->edge(temp_eh);
			if(edge.from_vertex().idx() == v1.idx() )
			{
				result[0].v1 = edge.to_vertex();
			}
			else if( edge.to_vertex().idx() == v1.idx() )
			{
				result[0].v1 = edge.from_vertex();
			}

			if(edge.from_vertex().idx() == v2.idx() )
			{
				result[1].v1 = edge.to_vertex();
			}
			else if( edge.to_vertex().idx() == v2.idx() )
			{
				result[1].v1 = edge.from_vertex();
			}
		}
	}

	face = mesh_->face(mesh_->face_handle(hf2));
	heh_Vec = face.halfedges();
	for( unsigned int i=0; i<heh_Vec.size(); ++i )
	{
		temp_eh = mesh_->edge_handle(heh_Vec[i]);
		if( temp_eh.idx() != eh.idx() )
		{
			edge = mesh_->edge(temp_eh);
			if(edge.from_vertex().idx() == v1.idx() )
			{
				result[1].v2 = edge.to_vertex();
			}
			else if( edge.to_vertex().idx() == v1.idx() )
			{
				result[1].v2 = edge.from_vertex();
			}

			if(edge.from_vertex().idx() == v2.idx() )
			{
				result[0].v2 = edge.to_vertex();
			}
			else if( edge.to_vertex().idx() == v2.idx() )
			{
				result[0].v2 = edge.from_vertex();
			}
		}
	}

	return result;
}

OpenVolumeMesh::HalfFaceHandle adjHFinCell(OpenVolumeMesh::HalfFaceHandle hf, OpenVolumeMesh::EdgeHandle eh, VolumeMesh* mesh_)
{
	OpenVolumeMesh::CellHandle ch = mesh_->incident_cell(hf);
	assert(ch.idx() >= 0);
	OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hff_Vec = cell.halffaces();
	OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hff_Vec[0]);
	OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge(eh);
	int vcount = 0;
	OpenVolumeMesh::HalfFaceHandle result;
	for(unsigned int i=0;i<hff_Vec.size();++i)
	{
		if( hff_Vec[i].idx() != hf.idx() )
		{
			hfv_it = mesh_->hfv_iter(hff_Vec[i]);
			vcount = 0;
			for(hfv_it;hfv_it;++hfv_it)
			{
				if(hfv_it->idx() == edge.from_vertex().idx() || hfv_it->idx() == edge.to_vertex().idx())
				{
					++vcount;
				}
			}
			if(vcount==2)
			{
				result = hff_Vec[i];
				break;
			}
		}
	}
	return result;
}

OpenVolumeMesh::EdgeHandle arbitraryOneAdjVerticalEHinCell(OpenVolumeMesh::HalfFaceHandle hf, VolumeMesh* mesh_)
{
	OpenVolumeMesh::CellHandle ch = mesh_->incident_cell(hf);
	OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hff_Vec = cell.halffaces();
	OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face( mesh_->face_handle(hf) );
	std::vector<OpenVolumeMesh::HalfEdgeHandle> heh_Vec = face.halfedges();
	std::vector<OpenVolumeMesh::EdgeHandle> edge_id(4);
	for(unsigned int i=0;i<heh_Vec.size();++i)
	{
		edge_id[i] = mesh_->edge_handle(heh_Vec[i]);
	}

	OpenVolumeMesh::EdgeHandle temp_eh;
	OpenVolumeMesh::EdgeHandle result;
	std::vector<bool> HFVisited(hff_Vec.size(),false);
	int v1[2] = {0,0};
	for(unsigned int j=0;j<edge_id.size();++j)
	{
		for(unsigned int i=0;i<hff_Vec.size();++i)
		{
			if( hff_Vec[i].idx() != hf.idx())
			{
				face = mesh_->face( mesh_->face_handle(hff_Vec[i]) );
				heh_Vec = face.halfedges();
				for( unsigned int k=0; k<heh_Vec.size(); ++k )
				{
					temp_eh = mesh_->edge_handle(heh_Vec[k]);
					if( temp_eh.idx() == edge_id[j].idx() )
					{
						OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge(temp_eh);
						v1[0] = edge.from_vertex().idx(); v1[1] = edge.to_vertex().idx();
						for(unsigned int l=0;l<heh_Vec.size();++l)
						{
							if( temp_eh.idx() != mesh_->edge_handle(heh_Vec[l]).idx() )
							{
								edge = mesh_->edge(mesh_->edge_handle(heh_Vec[l]));
								if(v1[0] == edge.from_vertex().idx() || v1[1] == edge.from_vertex().idx() 
									|| v1[0] == edge.to_vertex().idx() || v1[1] == edge.to_vertex().idx())
								{
									result = mesh_->edge_handle(heh_Vec[l]);
									return result;
								}
							}
						}
					}
				}
			}
		}
	}

	return result;
}

OpenVolumeMesh::EdgeHandle OneAdjVerticalEHinCell(OpenVolumeMesh::CellHandle ch, OpenVolumeMesh::HalfFaceHandle hf, OpenVolumeMesh::VertexHandle vh, VolumeMesh* mesh_)
{
	OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hff_Vec = cell.halffaces();
	OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face( mesh_->face_handle(hf) );
	std::vector<OpenVolumeMesh::HalfEdgeHandle> heh_Vec = face.halfedges();
	OpenVolumeMesh::EdgeHandle temp_eh;
	for(unsigned int i=0;i<heh_Vec.size();++i)
	{
		temp_eh = mesh_->edge_handle(heh_Vec[i]);
		OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge(temp_eh);
		if( edge.from_vertex().idx() == vh.idx() || edge.to_vertex().idx() == vh.idx() )
		{
			break;
		}
	}
	for(unsigned int i=0;i<hff_Vec.size();++i)
	{
		if( hff_Vec[i] != hf )
		{
			face = mesh_->face( mesh_->face_handle(hff_Vec[i]) );
			heh_Vec = face.halfedges();
			bool edge_OK = false;
			for(unsigned int j=0; j<heh_Vec.size(); ++j)
			{
				if( mesh_->edge_handle(heh_Vec[j]) == temp_eh )
				{
					edge_OK  =true;
					break;
				}
			}
			if(edge_OK)
			{
				break;
			}
		}
	}
	OpenVolumeMesh::EdgeHandle result;
	for(unsigned int j=0; j<heh_Vec.size(); ++j)
	{
		if( mesh_->edge_handle(heh_Vec[j]) != temp_eh )
		{
			OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge(mesh_->edge_handle(heh_Vec[j]));
			if( edge.from_vertex().idx() == vh.idx() || edge.to_vertex() == vh.idx() )
			{
				result = mesh_->edge_handle(heh_Vec[j]);
				break;
			}
		}
	}

	return result;
}

OpenVolumeMesh::HalfFaceHandle halffaceOncCellWithVertex(OpenVolumeMesh::CellHandle ch, OpenVolumeMesh::EdgeHandle eh, OpenVolumeMesh::VertexHandle vh, VolumeMesh* mesh_)
{
	OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hff_Vec = cell.halffaces();
	OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hff_Vec[0]);
	OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face( mesh_->face_handle(hff_Vec[0]) );
	std::vector<OpenVolumeMesh::HalfEdgeHandle> hef_Vec = face.halfedges();
	bool face_OK = false; bool vertex_OK = false; bool edge_OK = true;
	OpenVolumeMesh::HalfFaceHandle result;
	for(unsigned int i=0;i<hff_Vec.size();++i)
	{
		hfv_it = mesh_->hfv_iter(hff_Vec[i]);
		face = mesh_->face( mesh_->face_handle(hff_Vec[i]) );
		hef_Vec = face.halfedges();
		vertex_OK = false; edge_OK = true;
		for(unsigned int j=0;j<hef_Vec.size();++j,++hfv_it)
		{
			if( hfv_it->idx() == vh.idx() )
			{
				vertex_OK = true;
			}
			if( mesh_->edge_handle(hef_Vec[j]).idx() == eh.idx() )
			{
				edge_OK = false;
				break;
			}
		}

		if(vertex_OK && edge_OK)
		{
			result = hff_Vec[i];
			break;
		}
	}
	return result;
}

std::vector< VertexEdgeOnDoublet > adjEHInorderWithVertex(std::vector<OpenVolumeMesh::HalfFaceHandle>& hf_Vec, OpenVolumeMesh::VertexHandle vh, VolumeMesh* mesh_)
{
	OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter( hf_Vec[0] );
	std::vector<OpenVolumeMesh::VertexHandle> vh_Vec; int index = 0;
	for( hfv_it; hfv_it; ++hfv_it )
	{
		vh_Vec.push_back(*hfv_it);
		if( hfv_it->idx() == vh.idx() )
		{
			index = vh_Vec.size() - 1;
		}
	}

	OpenVolumeMesh::VertexHandle v1 = vh_Vec[ (index - 1 + vh_Vec.size())%vh_Vec.size() ];
	OpenVolumeMesh::VertexHandle v2 = vh_Vec[ (index + 1 + vh_Vec.size())%vh_Vec.size() ];
	OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face(mesh_->face_handle(hf_Vec[0]));
	std::vector<OpenVolumeMesh::HalfEdgeHandle> heh_Vec = face.halfedges();
	OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge(mesh_->edge_handle(heh_Vec[0]));
	std::vector< VertexEdgeOnDoublet > result(2);
	for( unsigned int i=0; i<heh_Vec.size(); ++i )
	{
		edge = mesh_->edge(mesh_->edge_handle(heh_Vec[i]));
		if( (edge.from_vertex() == vh && edge.to_vertex() == v1) || (edge.to_vertex() == vh && edge.from_vertex() == v1))
		{
			result[0]= VertexEdgeOnDoublet( v1, mesh_->edge_handle(heh_Vec[i]) );
		}
		else if((edge.from_vertex() == vh && edge.to_vertex() == v2) || (edge.to_vertex() == vh && edge.from_vertex() == v2))
		{
			result[1]= VertexEdgeOnDoublet( v2, mesh_->edge_handle(heh_Vec[i]) );
		}
	}
	std::vector<int> face_order(hf_Vec.size());
	face_order[0] = 0;

	std::vector<bool> visitedHF(hf_Vec.size(),false);
	visitedHF[0] = true;

	while(result.size() <hf_Vec.size() )
	{
		for( unsigned int i = 0; i< hf_Vec.size(); ++i )
		{
			if( !visitedHF[i] )
			{
				hfv_it = mesh_->hfv_iter( hf_Vec[i] );
				vh_Vec.clear();
				for( hfv_it; hfv_it; ++hfv_it )
				{
					vh_Vec.push_back(*hfv_it);
					if( hfv_it->idx() == vh.idx() )
					{
						index = vh_Vec.size() - 1;
					}
				}

				if( vh_Vec[(index - 1 + vh_Vec.size())%vh_Vec.size()] == v2 )
				{
					v1 = v2; v2 = vh_Vec[(index + 1 + vh_Vec.size())%vh_Vec.size()];
					face = mesh_->face(mesh_->face_handle(hf_Vec[i]));
					heh_Vec = face.halfedges();
					bool add_edge_ok = false;
					for( unsigned int j=0; j<heh_Vec.size(); ++j )
					{
						edge = mesh_->edge(mesh_->edge_handle(heh_Vec[j]));
						if((edge.from_vertex() == vh && edge.to_vertex() == v2) || (edge.to_vertex() == vh && edge.from_vertex() == v2))
						{
							face_order[ result.size() -1 ] = i;
							result.push_back( VertexEdgeOnDoublet( v2, mesh_->edge_handle(heh_Vec[j]) ) );
							add_edge_ok = true;
							visitedHF[i] = true;
							break;
						}
					}
					if(add_edge_ok)
					{break;}
				}
			}
		}
	}

	for(unsigned int i=0;i<visitedHF.size();++i)
	{
		if(!visitedHF[i])
		{
			face_order[result.size() - 1] = i;
			break;
		}
	}

	std::vector<OpenVolumeMesh::HalfFaceHandle> temp_hf;
	for(unsigned int i=0;i<face_order.size();++i)
	{
		temp_hf.push_back(hf_Vec[face_order[i]]);
	}

	hf_Vec = temp_hf;

	return result;

}

double compute_average_edge_len(VolumeMesh* mesh_)
{
	OpenVolumeMesh::EdgeIter e_it = mesh_->edges_begin();
	double all_len = 0.0;
	OpenVolumeMesh::Geometry::Vec3d p_from;
	OpenVolumeMesh::Geometry::Vec3d p_to;
	for(e_it;e_it != mesh_->edges_end();++e_it)
	{
		OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge(*e_it);
		p_from = mesh_->vertex(edge.from_vertex());
		p_to   = mesh_->vertex(edge.to_vertex());
		all_len += (p_from - p_to).norm();
	}
	return all_len/mesh_->n_edges();
}

double compute_average_face_area(VolumeMesh* mesh_)//for tet mesh
{
	OpenVolumeMesh::FaceIter f_it = mesh_->faces_begin();
	double all_area = 0.0;
	OpenVolumeMesh::Geometry::Vec3d p_0;
	OpenVolumeMesh::Geometry::Vec3d p_1;
	OpenVolumeMesh::Geometry::Vec3d p_2;
	double count_face =0 ;
	for(f_it;f_it != mesh_->faces_end();++f_it)
	{
		OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(mesh_->halfface_handle(*f_it,0));
		p_0 = mesh_->vertex(*hfv_it); ++hfv_it;
		p_1 = mesh_->vertex(*hfv_it); ++hfv_it;
		p_2 = mesh_->vertex(*hfv_it);
		
		all_area += (OpenVolumeMesh::Geometry::cross(p_0 - p_1, p_2 - p_1)).norm();
		count_face += 1.0;
	}
	return (all_area*0.5/count_face);
}

double compute_absolute_volume(VolumeMesh* mesh_, OpenVolumeMesh::CellHandle cell_handle)
{
	Eigen::Matrix4d volume_matrix;
	OpenVolumeMesh::CellVertexIter cv_it = mesh_->cv_iter(cell_handle);
	OpenVolumeMesh::Geometry::Vec3d P; int count = 0;
	for(cv_it;cv_it.valid();++cv_it)
	{
		P = mesh_->vertex(*cv_it);
		volume_matrix(1,count) = P[0];
		volume_matrix(2,count) = P[1];
		volume_matrix(3,count) = P[2];
		volume_matrix(0,count) = 1.0;
	}

	return std::abs( volume_matrix.determinant() ) / 6.0;
}

void compute_all_absolute_volume(VolumeMesh* mesh_, std::vector<double>& cell_vol)
{
	cell_vol.resize(mesh_->n_cells());
	OpenVolumeMesh::CellIter c_it = mesh_->cells_begin();
	OpenVolumeMesh::CellVertexIter cv_it = mesh_->cv_iter(*c_it);
	Eigen::Matrix4d volume_matrix;
	OpenVolumeMesh::Geometry::Vec3d P; int count = 0;
	for(c_it;c_it != mesh_->cells_end();++c_it)
	{
		cv_it = mesh_->cv_iter(*c_it);
		count = 0;
		for( cv_it; cv_it.valid(); ++cv_it )
		{
			P = mesh_->vertex(*cv_it);
			volume_matrix(1,count) = P[0];
			volume_matrix(2,count) = P[1];
			volume_matrix(3,count) = P[2];
			volume_matrix(0,count) = 1.0;
			++count;
		}

		cell_vol[c_it->idx()] = std::abs( volume_matrix.determinant() ) / 6.0;
	}
}

//false : have negative volume
bool check_negative_volume(VolumeMesh* mesh_)
{
	unsigned nc = mesh_->n_cells(); std::vector<int> flag_cell(nc, -1);
	OpenVolumeMesh::Geometry::Vec3d p, p1, p2, p3; double min_volume = 1e30; int min_cell_id = -1;
	for(OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
	{
		int v_id = v_it->idx(); p = mesh_->vertex(*v_it);
		for(OpenVolumeMesh::VertexCellIter vc_it = mesh_->vc_iter(*v_it); vc_it; ++vc_it)
		{
			int c_id = vc_it->idx();
			if(flag_cell[c_id] == -1)
			{
				OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell( vc_it.cur_handle() );
				std::vector<OpenVolumeMesh::HalfFaceHandle> hfh_vec = cell.halffaces();
				for(unsigned i = 0; i < hfh_vec.size(); ++i)
				{
					bool find_hf = true;
					for(OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfh_vec[i]); hfv_it; ++hfv_it)
					{
						if(hfv_it->idx() == v_id)
						{
							find_hf = false;
							break;
						}
					}
					if( find_hf )
					{
						OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfh_vec[i]); 
						p1 = mesh_->vertex( hfv_it.cur_handle() );
						++hfv_it; p2 = mesh_->vertex( hfv_it.cur_handle() );
						++hfv_it; p3 = mesh_->vertex( hfv_it.cur_handle() );
						break;
					}
				}

				double c_v = OpenVolumeMesh::Geometry::dot(p - p1, OpenVolumeMesh::Geometry::cross(p2 - p1, p3 - p1));
				if(c_v < min_volume) {min_volume = c_v; min_cell_id = c_id;}

				flag_cell[c_id] = 1;
			}
		}
	}
	printf("Min Volume : %20.19f on %d\n", min_volume, min_cell_id);
	if(min_volume < 1e-16)
	{
		return false;
	}
	return true;
}

//false : have negative volume
bool check_negative_volume_id(VolumeMesh* mesh_, int& min_cell_id)
{
	unsigned nc = mesh_->n_cells(); std::vector<int> flag_cell(nc, -1);
	OpenVolumeMesh::Geometry::Vec3d p, p1, p2, p3; double min_volume = 1e30;
	for(OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
	{
		int v_id = v_it->idx(); p = mesh_->vertex(*v_it);
		for(OpenVolumeMesh::VertexCellIter vc_it = mesh_->vc_iter(*v_it); vc_it; ++vc_it)
		{
			int c_id = vc_it->idx();
			if(flag_cell[c_id] == -1)
			{
				OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell( vc_it.cur_handle() );
				std::vector<OpenVolumeMesh::HalfFaceHandle> hfh_vec = cell.halffaces();
				for(unsigned i = 0; i < hfh_vec.size(); ++i)
				{
					bool find_hf = true;
					for(OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfh_vec[i]); hfv_it; ++hfv_it)
					{
						if(hfv_it->idx() == v_id)
						{
							find_hf = false;
							break;
						}
					}
					if( find_hf )
					{
						OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfh_vec[i]); 
						p1 = mesh_->vertex( hfv_it.cur_handle() );
						++hfv_it; p2 = mesh_->vertex( hfv_it.cur_handle() );
						++hfv_it; p3 = mesh_->vertex( hfv_it.cur_handle() );
						break;
					}
				}

				double c_v = OpenVolumeMesh::Geometry::dot(p - p1, OpenVolumeMesh::Geometry::cross(p2 - p1, p3 - p1));
				if(c_v < min_volume) {min_volume = c_v; min_cell_id = c_id;}

				flag_cell[c_id] = 1;
			}
		}
	}
	printf("Min Volume : %20.19f on %d\n", min_volume, min_cell_id);
	if(min_volume < 1e-16)
	{
		return false;
	}
	return true;
}

void get_negative_volume(VolumeMesh* mesh_, std::vector<OpenVolumeMesh::CellHandle>& negative_cell)
{
	negative_cell.clear();
	unsigned nc = mesh_->n_cells(); std::vector<int> flag_cell(nc, -1);
	OpenVolumeMesh::Geometry::Vec3d p, p1, p2, p3; double min_volume = 1e30; int min_cell_id = -1;
	for(OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
	{
		int v_id = v_it->idx(); p = mesh_->vertex(*v_it);
		for(OpenVolumeMesh::VertexCellIter vc_it = mesh_->vc_iter(*v_it); vc_it; ++vc_it)
		{
			int c_id = vc_it->idx();
			if(flag_cell[c_id] == -1)
			{
				OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell( vc_it.cur_handle() );
				std::vector<OpenVolumeMesh::HalfFaceHandle> hfh_vec = cell.halffaces();
				for(unsigned i = 0; i < hfh_vec.size(); ++i)
				{
					bool find_hf = true;
					for(OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfh_vec[i]); hfv_it; ++hfv_it)
					{
						if(hfv_it->idx() == v_id)
						{
							find_hf = false;
							break;
						}
					}
					if( find_hf )
					{
						OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfh_vec[i]); 
						p1 = mesh_->vertex( hfv_it.cur_handle() );
						++hfv_it; p2 = mesh_->vertex( hfv_it.cur_handle() );
						++hfv_it; p3 = mesh_->vertex( hfv_it.cur_handle() );
						break;
					}
				}

				double c_v = OpenVolumeMesh::Geometry::dot(p - p1, OpenVolumeMesh::Geometry::cross(p2 - p1, p3 - p1));
				if(c_v < 0.0) { negative_cell.push_back( vc_it.cur_handle() ); }
				flag_cell[c_id] = 1;
			}
		}
	}
	//printf("Negative Cell : %d\n", negative_cell.size());
}

void calc_tet_face_normal(VolumeMesh* mesh_, OpenVolumeMesh::HalfFaceHandle& hf, OpenVolumeMesh::Geometry::Vec3d& n)
{
	OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hf);
	OpenVolumeMesh::Geometry::Vec3d p0 = mesh_->vertex(*hfv_it);
	++hfv_it; OpenVolumeMesh::Geometry::Vec3d p1 = mesh_->vertex(*hfv_it);
	++hfv_it; OpenVolumeMesh::Geometry::Vec3d p2 = mesh_->vertex(*hfv_it);
	n = ( OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0 ) ).normalize();
}

void calc_tet_face_normal_anisotropy(VolumeMesh* mesh_, OpenVolumeMesh::HalfFaceHandle& hf,const OpenVolumeMesh::Geometry::Vec6d& Q, OpenVolumeMesh::Geometry::Vec3d& n)
{
	OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hf);
	OpenVolumeMesh::Geometry::Vec3d p0_ = mesh_->vertex(*hfv_it);
	++hfv_it; OpenVolumeMesh::Geometry::Vec3d p1_ = mesh_->vertex(*hfv_it);
	++hfv_it; OpenVolumeMesh::Geometry::Vec3d p2_ = mesh_->vertex(*hfv_it);
	OpenVolumeMesh::Geometry::Vec3d p0(Q[0]*p0_[0] + Q[1]*p0_[1] + Q[2]*p0_[2], Q[1]*p0_[0] + Q[3]*p0_[1] + Q[4]*p0_[2], Q[2]*p0_[0] + Q[4]*p0_[1] + Q[5]*p0_[2]);
	OpenVolumeMesh::Geometry::Vec3d p1(Q[0]*p1_[0] + Q[1]*p1_[1] + Q[2]*p1_[2], Q[1]*p1_[0] + Q[3]*p1_[1] + Q[4]*p1_[2], Q[2]*p1_[0] + Q[4]*p1_[1] + Q[5]*p1_[2]);
	OpenVolumeMesh::Geometry::Vec3d p2(Q[0]*p2_[0] + Q[1]*p2_[1] + Q[2]*p2_[2], Q[1]*p2_[0] + Q[3]*p2_[1] + Q[4]*p2_[2], Q[2]*p2_[0] + Q[4]*p2_[1] + Q[5]*p2_[2]);
	n = ( OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0 ) ).normalize();
}

//false : not in the tet
bool check_in_tet(VolumeMesh* mesh_, const OpenVolumeMesh::Geometry::Vec3d& p, OpenVolumeMesh::CellHandle ch)
{
	OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(ch);
	std::vector<OpenVolumeMesh::HalfFaceHandle>& hfh_vec = cell.get_halffaces();
	std::vector<OpenVolumeMesh::Geometry::Vec3d> one_face(3); int count = 0;
	for (int i = 0; i < hfh_vec.size(); ++i)
	{
		count = 0;
		for (OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfh_vec[i]); hfv_it; ++hfv_it)
		{
			one_face[count] = mesh_->vertex(*hfv_it); ++count;
		}
		double v = OpenVolumeMesh::Geometry::dot(p - one_face[0], OpenVolumeMesh::Geometry::cross(one_face[1] - one_face[0], one_face[2] - one_face[0]));
		if (v < 0)
		{
			return false;
		}
	}
	return true;
}

void compute_point_area(SurfaceMesh* mesh_, std::vector<std::map<int,double>>& cornerArea, std::vector<double>& pointArea)
{
	pointArea.resize(mesh_->n_vertices(),0.0);
	cornerArea.resize(mesh_->n_faces());

	SurfaceMesh::FaceIter f_it = mesh_->faces_begin();
	int temp_face_index = 0;
	SurfaceMesh::FaceHalfedgeIter fhe_it;
	SurfaceMesh::Point e[3];
	int v[3];

	//#pragma omp parallel for
	for( f_it; f_it != mesh_->faces_end(); ++f_it)
	{
		temp_face_index = f_it.handle().idx();

		fhe_it = mesh_->fh_iter(f_it);
		e[0] = mesh_->point( mesh_->to_vertex_handle(fhe_it) ) - mesh_->point( mesh_->from_vertex_handle(fhe_it) );
		v[2] = mesh_->to_vertex_handle(fhe_it).idx();
		++fhe_it;
		e[1] = mesh_->point( mesh_->to_vertex_handle(fhe_it) ) - mesh_->point( mesh_->from_vertex_handle(fhe_it) );
		v[0] = mesh_->to_vertex_handle(fhe_it).idx();
		++fhe_it;
		e[2] = mesh_->point( mesh_->to_vertex_handle(fhe_it) ) - mesh_->point( mesh_->from_vertex_handle(fhe_it) );
		v[1] = mesh_->to_vertex_handle(fhe_it).idx();

		double area = 0.5f * ( e[0]%e[1] ).norm();
		double l2[3] = { e[0].sqrnorm(), e[1].sqrnorm(), e[2].sqrnorm() };
		double ew[3] = { l2[0] * (l2[1] + l2[2] - l2[0]),
			l2[1] * (l2[2] + l2[0] - l2[1]),
			l2[2] * (l2[0] + l2[1] - l2[2]) };

		if (ew[0] <= 0.0) 
		{
			cornerArea[temp_face_index][v[1]] = -0.25 * l2[2] * area / OpenMesh::dot(e[0],e[2]);
			cornerArea[temp_face_index][v[2]] = -0.25 * l2[1] * area / OpenMesh::dot(e[0],e[1]);
			cornerArea[temp_face_index][v[0]] = area - cornerArea[temp_face_index][v[1]] -cornerArea[temp_face_index][v[2]];
		}
		else if (ew[1] <= 0.0)
		{
			cornerArea[temp_face_index][v[2]] = -0.25 * l2[0] * area / OpenMesh::dot(e[1],e[0]);
			cornerArea[temp_face_index][v[0]] = -0.25 * l2[2] * area / OpenMesh::dot(e[1],e[2]);
			cornerArea[temp_face_index][v[1]] = area - cornerArea[temp_face_index][v[2]] - cornerArea[temp_face_index][v[0]];
		}
		else if (ew[2] <= 0.0f)
		{
			cornerArea[temp_face_index][v[0]] = -0.25 * l2[1] * area / OpenMesh::dot(e[2],e[1]);
			cornerArea[temp_face_index][v[1]] = -0.25 * l2[0] * area / OpenMesh::dot(e[2],e[0]);
			cornerArea[temp_face_index][v[2]] = area - cornerArea[temp_face_index][v[0]] - cornerArea[temp_face_index][v[1]];
		}
		else 
		{
			double ewscale = 0.5 * area / (ew[0] + ew[1] + ew[2]);
			for (int j = 0; j < 3; j++)
			{
				cornerArea[temp_face_index][v[j]] = ewscale * (ew[(j+1)%3] + ew[(j+2)%3]);
			}
		}

		//#pragma omp atomic
		pointArea[v[0]] += cornerArea[temp_face_index][v[0]];
		//#pragma omp atomic
		pointArea[v[1]] += cornerArea[temp_face_index][v[1]];
		//#pragma omp atomic
		pointArea[v[2]] += cornerArea[temp_face_index][v[2]];
	}
}

void rot_coord_sys(const OpenMesh::Vec3d &old_u, const OpenMesh::Vec3d &old_v, const OpenMesh::Vec3d &new_norm, OpenMesh::Vec3d &new_u, OpenMesh::Vec3d &new_v)
{
	new_u = old_u;
	new_v = old_v;
	OpenMesh::Vec3d old_norm = OpenMesh::cross( old_u , old_v );
	double ndot = OpenMesh::dot( old_norm , new_norm);
	if ( ndot <= -1.0 )
	{
		new_u = -new_u;
		new_v = -new_v;
		return;
	}
	OpenMesh::Vec3d perp_old = new_norm - ndot * old_norm;
	OpenMesh::Vec3d dperp = 1.0f / (1 + ndot) * (old_norm + new_norm);
	new_u -= dperp * OpenMesh::dot(new_u , perp_old);
	new_v -= dperp * OpenMesh::dot(new_v , perp_old);
}

void proj_curv(const OpenMesh::Vec3d &old_u, const OpenMesh::Vec3d &old_v, double old_ku, double old_kuv, double old_kv, const OpenMesh::Vec3d &new_u, const OpenMesh::Vec3d &new_v, double &new_ku, double &new_kuv, double &new_kv)
{
	OpenMesh::Vec3d r_new_u; OpenMesh::Vec3d r_new_v;
	rot_coord_sys(new_u, new_v, OpenMesh::cross(old_u, old_v), r_new_u, r_new_v);

	double u1 = OpenMesh::dot( r_new_u , old_u);
	double v1 = OpenMesh::dot( r_new_u , old_v);
	double u2 = OpenMesh::dot( r_new_v , old_u);
	double v2 = OpenMesh::dot( r_new_v , old_v);
	new_ku  = old_ku * u1*u1 + old_kuv * (2.0f  * u1*v1) + old_kv * v1*v1;
	new_kuv = old_ku * u1*u2 + old_kuv * (u1*v2 + u2*v1) + old_kv * v1*v2;
	new_kv  = old_ku * u2*u2 + old_kuv * (2.0f  * u2*v2) + old_kv * v2*v2;
}

void diagonalize_curv(const OpenMesh::Vec3d &old_u, const OpenMesh::Vec3d &old_v, double ku, double kuv, double kv, const OpenMesh::Vec3d &new_norm, OpenMesh::Vec3d &pdir1, OpenMesh::Vec3d &pdir2, double &vk1, double &vk2)
{
	OpenMesh::Vec3d r_old_u, r_old_v;
	rot_coord_sys(old_u, old_v, new_norm, r_old_u, r_old_v);

	double c = 1.0, s = 0.0, tt = 0.0;
	if(kuv != 0.0) 
	{
		// Jacobi rotation to diagonalize
		double h = 0.5 * (kv - ku) / kuv;
		tt = (h < 0.0) ?
			1.0 / (h - std::sqrt(1.0 + h*h)) :
		1.0 / (h + std::sqrt(1.0 + h*h));
		c = 1.0 / std::sqrt(1.0 + tt*tt);
		s = tt * c;
	}

	vk1 = ku - tt * kuv;
	vk2 = kv + tt * kuv;

	if (std::abs(vk1) >= std::abs(vk2))
	{
		pdir1 = c*r_old_u - s*r_old_v;
	} 
	else
	{
		std::swap(vk1, vk2);
		pdir1 = s*r_old_u + c*r_old_v;
	}
	pdir2 = OpenMesh::cross( new_norm , pdir1);
}

void compute_principal_curvature(SurfaceMesh* mesh_, 
								 std::vector<double>& K1, std::vector<double>& K2, 
								 std::vector<OpenMesh::Vec3d>& dir1,std::vector<OpenMesh::Vec3d>& dir2)
{
	if(!mesh_->has_vertex_normals())
	{
		mesh_->request_vertex_normals();
		mesh_->update_vertex_normals();
	}

	int nv = mesh_->n_vertices();
	int nf = mesh_->n_faces();
	SurfaceMesh::FaceIter f_it;

	//compute vertex normal
	std::vector<OpenMesh::Vec3d> vertex_normal(nv,OpenMesh::Vec3d(0.0,0.0,0.0));
	SurfaceMesh::FaceVertexIter fv_it;
	OpenMesh::Vec3d v[3];
	int v_index[3];
	for(f_it = mesh_->faces_begin();f_it != mesh_->faces_end(); ++f_it)
	{
		fv_it = mesh_->fv_iter(f_it);
		int count = 0;
		for(fv_it;fv_it;++fv_it)
		{
			v[count] = mesh_->point(fv_it);
			v_index[count] = fv_it.handle().idx();
			++count;
		}
		OpenMesh::Vec3d a = v[0] - v[1];
		OpenMesh::Vec3d b = v[1] - v[2];
		OpenMesh::Vec3d c = v[2] - v[0];
		double l2a = a.sqrnorm(); double l2b = b.sqrnorm(); double l2c = c.sqrnorm();
		if (!l2a || !l2b || !l2c)
			continue;

		OpenMesh::Vec3d temp_fn = OpenMesh::cross(a,b);
		vertex_normal[v_index[0]] += temp_fn * (1.0 / (l2a * l2c));
		vertex_normal[v_index[1]] += temp_fn * (1.0 / (l2b * l2a));
		vertex_normal[v_index[2]] += temp_fn * (1.0 / (l2c * l2b));
	}

	K1.clear(); K1.resize(nv,0.0); K2.clear(); K2.resize(nv,0.0);
	dir1.clear(); dir1.resize(nv,OpenMesh::Vec3d(0,0,0));
	dir2.clear(); dir2.resize(nv,OpenMesh::Vec3d(0,0,0));
	std::vector<double> k12(nv,0.0);

	SurfaceMesh::FaceHalfedgeIter fhe_it;
	for(f_it = mesh_->faces_begin();f_it != mesh_->faces_end(); ++f_it)
	{
		fhe_it = mesh_->fh_iter(f_it);
		for(fhe_it;fhe_it;++fhe_it)
		{
			dir1[mesh_->from_vertex_handle(fhe_it).idx()] =
				mesh_->point( mesh_->to_vertex_handle(fhe_it) ) 
				- mesh_->point( mesh_->from_vertex_handle(fhe_it));
		}
	}
	
	SurfaceMesh::VertexIter v_it;
	int vertex_id; OpenMesh::Vec3d vn;
	for(v_it = mesh_->vertices_begin();v_it != mesh_->vertices_end(); ++v_it)
	{
		//vn = mesh_->normal(v_it); vn.normalize();
		vertex_normal[v_it.handle().idx()].normalize();
		vn = vertex_normal[v_it.handle().idx()];
		vertex_id = v_it.handle().idx();
		dir1[ vertex_id ] = OpenMesh::cross(dir1[ vertex_id ], vn);
		dir1[ vertex_id ].normalize();
		dir2[ vertex_id ] = OpenMesh::cross(vn, dir1[ vertex_id ]);
		dir2[ vertex_id ].normalize();
	}

	std::vector<std::map<int,double>> cornerArea;
	std::vector<double> pointArea;
	compute_point_area(mesh_, cornerArea, pointArea);

	OpenMesh::Vec3d ev[3]; 
	SurfaceMesh::HalfedgeHandle heh; SurfaceMesh::HalfedgeHandle temp_heh; 
	OpenMesh::Vec3d t; OpenMesh::Vec3d n;OpenMesh::Vec3d b;
	for(f_it = mesh_->faces_begin();f_it != mesh_->faces_end(); ++f_it)
	{
		// Edges
		fhe_it = mesh_->fh_iter(f_it); heh = fhe_it.handle();
		temp_heh = mesh_->next_halfedge_handle(heh);
		ev[0] = mesh_->point(mesh_->to_vertex_handle(temp_heh)) - mesh_->point(mesh_->from_vertex_handle(temp_heh));
		temp_heh = mesh_->prev_halfedge_handle(heh);
		ev[1] = mesh_->point(mesh_->to_vertex_handle(temp_heh)) - mesh_->point(mesh_->from_vertex_handle(temp_heh));
		ev[2] = mesh_->point(mesh_->to_vertex_handle(     heh)) - mesh_->point(mesh_->from_vertex_handle(     heh));
		
		// N-T-B coordinate system per face
		t = ev[0]; t.normalize();
		n = OpenMesh::cross(ev[0],ev[1]);
		b = OpenMesh::cross(n,t); b.normalize();

		// Estimate curvature based on variation of normals
		// along edges
		temp_heh = mesh_->next_halfedge_handle(heh);
		Eigen::Vector3d m(0.0,0.0,0.0);
		Eigen::Vector3d x;
		Eigen::Matrix3d w; w.setZero();
		for(int j=0;j<3;++j)
		{
			double u = OpenMesh::dot(ev[j],t);
			double v = OpenMesh::dot(ev[j],b);
			w(0,0) += u*u;
			w(0,1) += u*v;
			//w[1][1] += v*v + u*u; 
			//w[1][2] += u*v; 
			w(2,2) += v*v;
			/*OpenMesh::Vec3d dn = mesh_->normal(mesh_->to_vertex_handle(temp_heh)) 
			- mesh_->normal(mesh_->from_vertex_handle(temp_heh));*/
			OpenMesh::Vec3d dn = vertex_normal[mesh_->to_vertex_handle(temp_heh).idx()] 
							   - vertex_normal[mesh_->from_vertex_handle(temp_heh).idx()];
			double dnu = OpenMesh::dot(dn, t);
			double dnv = OpenMesh::dot(dn, b);
			m(0) += dnu*u;
			m(1) += dnu*v + dnv*u;
			m(2) += dnv*v;
			temp_heh = mesh_->next_halfedge_handle(temp_heh);
		}
		w(1,1) = w(0,0) + w(2,2);
		w(1,2) = w(0,1);
		w(2,1) = w(1,2);
		w(1,0) = w(0,1);
		//std::cout << w;
		x = w.fullPivHouseholderQr().solve(m);

		temp_heh = heh;
		for(int j=0;j<3;++j)
		{
			vertex_id = mesh_->from_vertex_handle(temp_heh).idx();
			double c1, c12, c2;
			proj_curv(t, b, x(0), x(1), x(2),
				dir1[vertex_id], dir2[vertex_id], c1, c12, c2);
			double wt = cornerArea[f_it.handle().idx()][vertex_id] / pointArea[vertex_id];
			K1[vertex_id]  += wt * c1;
			k12[vertex_id] += wt * c12;
			K2[vertex_id] += wt * c2;
			temp_heh = mesh_->next_halfedge_handle(temp_heh);
		}
	}

	for(v_it = mesh_->vertices_begin();v_it != mesh_->vertices_end(); ++v_it)
	{
		vertex_id = v_it.handle().idx();
		/*diagonalize_curv(dir1[vertex_id], dir2[vertex_id],
		k1[vertex_id], k12[vertex_id], k2[vertex_id],
		mesh_->normal(v_it), dir1[vertex_id], dir2[vertex_id],
		k1[vertex_id], k2[vertex_id]);*/
		diagonalize_curv(dir1[vertex_id], dir2[vertex_id],
			K1[vertex_id], k12[vertex_id], K2[vertex_id],
			vertex_normal[vertex_id], dir1[vertex_id], dir2[vertex_id],
			K1[vertex_id], K2[vertex_id]);
	}
}

bool line_intersect_triangle(const OpenVolumeMesh::Geometry::Vec3d& p,
							 const OpenVolumeMesh::Geometry::Vec3d& q, 
							 const OpenVolumeMesh::Geometry::Vec3d& r, 
							 const OpenVolumeMesh::Geometry::Vec3d& s, 
							 const OpenVolumeMesh::Geometry::Vec3d& t)
{
	OpenVolumeMesh::Geometry::Vec3d Normal = OpenVolumeMesh::Geometry::cross(s-r, t-r);
	double s0 = OpenVolumeMesh::Geometry::dot( Normal, p - r );
	double s1 = OpenVolumeMesh::Geometry::dot( Normal, q - r );
	if (s0 * s1 >= 0)
	{
		return false;
	}

	double temp = -s0 / (s1 - s0);

	if (temp <= 0 || temp >= 1)
	{
		return false;
	}

	OpenVolumeMesh::Geometry::Vec3d O = p + temp * (q - p);
	OpenVolumeMesh::Geometry::Vec3d C[3];
	C[0] = O - r, C[1] = O - s, C[2] = O - t;
	double u = OpenVolumeMesh::Geometry::dot(Normal, OpenVolumeMesh::Geometry::cross(C[0], C[1]));
	double v = OpenVolumeMesh::Geometry::dot(Normal, OpenVolumeMesh::Geometry::cross(C[1], C[2]));
	double w = OpenVolumeMesh::Geometry::dot(Normal, OpenVolumeMesh::Geometry::cross(C[2], C[0]));
	double uvw = u + v + w;
	u /= uvw, v /= uvw;

	if (u < 1.0e-8 || v < 1.0e-8 || u + v > 1 - 1.0e-8)
		//if (u < -1.0e-8 || v < -1.0e-8 || u + v > 1 + 1.0e-8)
	{
		return false;
	}

	return true;
}

bool is_flip_face_23_ok(VolumeMesh* mesh_, OpenVolumeMesh::FaceHandle fh)
{
	OpenVolumeMesh::HalfFaceHandle hfh0 = mesh_->halfface_handle(fh, 0);
	OpenVolumeMesh::HalfFaceHandle hfh1 = mesh_->halfface_handle(fh, 1);
	
	OpenVolumeMesh::CellHandle ch0 = mesh_->incident_cell(hfh0);
	OpenVolumeMesh::CellHandle ch1 = mesh_->incident_cell(hfh1);

	if(ch0 == VolumeMesh::InvalidCellHandle || ch1 == VolumeMesh::InvalidCellHandle)
	{
		return false;
	}

	OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfh0);
	OpenVolumeMesh::VertexHandle vh = hfv_it.cur_handle();
	OpenVolumeMesh::Geometry::Vec3d r = mesh_->vertex(vh); int r_id = vh.idx();
	++hfv_it; vh = hfv_it.cur_handle();
	OpenVolumeMesh::Geometry::Vec3d s = mesh_->vertex(vh); int s_id = vh.idx();
	++hfv_it; vh = hfv_it.cur_handle();
	OpenVolumeMesh::Geometry::Vec3d t = mesh_->vertex(vh); int t_id = vh.idx();
	OpenVolumeMesh::Geometry::Vec3d p, q; OpenVolumeMesh::VertexHandle p_vh, q_vh;
	for(OpenVolumeMesh::CellVertexIter cv_it = mesh_->cv_iter(ch0); cv_it; ++cv_it)
	{
		int v_id = cv_it->idx();
		if( v_id != r_id && v_id != s_id && v_id != t_id )
		{
			p_vh = cv_it.cur_handle();
			p = mesh_->vertex( p_vh );
			break;
		}
	}
	for(OpenVolumeMesh::CellVertexIter cv_it = mesh_->cv_iter(ch1); cv_it; ++cv_it)
	{
		int v_id = cv_it->idx();
		if( v_id != r_id && v_id != s_id && v_id != t_id )
		{
			q_vh = cv_it.cur_handle();
			q = mesh_->vertex( q_vh );
			break;
		}
	}

	for(OpenVolumeMesh::VertexOHalfEdgeIter voh_it = mesh_->voh_iter(p_vh); voh_it; ++voh_it)
	{
		OpenVolumeMesh::OpenVolumeMeshEdge e = mesh_->edge(mesh_->edge_handle(voh_it.cur_handle()));
		if( (e.from_vertex() == p_vh && e.to_vertex() == q_vh) || (e.from_vertex() == q_vh && e.to_vertex() == p_vh) )
		{
			return false;
		}
	}

	if( !line_intersect_triangle(p, q, r, s, t) )
	{
		return false;
	}

	return true;
}

bool flip_face_23(VolumeMesh* mesh_, OpenVolumeMesh::FaceHandle fh)
{
	if( !is_flip_face_23_ok(mesh_, fh) )
	{return false;}

	OpenVolumeMesh::HalfFaceHandle hfh0 = mesh_->halfface_handle(fh, 0);
	OpenVolumeMesh::HalfFaceHandle hfh1 = mesh_->halfface_handle(fh, 1);

	OpenVolumeMesh::CellHandle ch0 = mesh_->incident_cell(hfh0);
	OpenVolumeMesh::CellHandle ch1 = mesh_->incident_cell(hfh1);

	if(ch0 == VolumeMesh::InvalidCellHandle || ch1 == VolumeMesh::InvalidCellHandle)
	{
		return false;
	}

	OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfh0);
	OpenVolumeMesh::VertexHandle vh = hfv_it.cur_handle();
	OpenVolumeMesh::Geometry::Vec3d r = mesh_->vertex(vh); int r_id = vh.idx();
	++hfv_it; vh = hfv_it.cur_handle();
	OpenVolumeMesh::Geometry::Vec3d s = mesh_->vertex(vh); int s_id = vh.idx();
	++hfv_it; vh = hfv_it.cur_handle();
	OpenVolumeMesh::Geometry::Vec3d t = mesh_->vertex(vh); int t_id = vh.idx();
	OpenVolumeMesh::Geometry::Vec3d p, q; OpenVolumeMesh::VertexHandle p_vh, q_vh; int p_id, q_id;
	for(OpenVolumeMesh::CellVertexIter cv_it = mesh_->cv_iter(ch0); cv_it; ++cv_it)
	{
		int v_id = cv_it->idx();
		if( v_id != r_id && v_id != s_id && v_id != t_id )
		{
			p_vh = cv_it.cur_handle(); p_id = p_vh.idx();
			p = mesh_->vertex( p_vh );
			break;
		}
	}
	for(OpenVolumeMesh::CellVertexIter cv_it = mesh_->cv_iter(ch1); cv_it; ++cv_it)
	{
		int v_id = cv_it->idx();
		if( v_id != r_id && v_id != s_id && v_id != t_id )
		{
			q_vh = cv_it.cur_handle(); q_id = q_vh.idx();
			q = mesh_->vertex( q_vh );
			break;
		}
	}

	OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(ch0);
	std::vector<OpenVolumeMesh::HalfFaceHandle> hfh_vec = cell.halffaces();
	OpenVolumeMesh::HalfFaceHandle hf_prt, hf_psr, hf_pts; bool prt_ok = false; bool psr_ok = false; bool pts_ok = false;
	for(unsigned i=0; i<hfh_vec.size(); ++i)
	{
		if( hfh_vec[i] != hfh0 )
		{
			OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfh_vec[i]); int v_id0 = hfv_it->idx();
			++hfv_it; int v_id1 = hfv_it->idx(); ++hfv_it; int v_id2 = hfv_it->idx();
			if( !prt_ok )
			{
				if( (v_id0 == p_id && v_id1 == r_id && v_id2 == t_id)
					|| (v_id1 == p_id && v_id2 == r_id && v_id0 == t_id)
					|| (v_id2 == p_id && v_id0 == r_id && v_id1 == t_id) )
				{
					hf_prt = hfh_vec[i];
					prt_ok = true;
				}
			}

			if( !psr_ok )
			{
				if( (v_id0 == p_id && v_id1 == s_id && v_id2 == r_id)
					|| (v_id1 == p_id && v_id2 == s_id && v_id0 == r_id)
					|| (v_id2 == p_id && v_id0 == s_id && v_id1 == r_id) )
				{
					hf_psr = hfh_vec[i];
					psr_ok = true;
				}
			}

			if( !pts_ok )
			{
				if( (v_id0 == p_id && v_id1 == t_id && v_id2 == s_id)
					|| (v_id1 == p_id && v_id2 == t_id && v_id0 == s_id)
					|| (v_id2 == p_id && v_id0 == t_id && v_id1 == s_id) )
				{
					hf_pts = hfh_vec[i];
					pts_ok = true;
				}
			}
		}
	}

	cell = mesh_->cell(ch1);
	hfh_vec = cell.halffaces();
	OpenVolumeMesh::HalfFaceHandle hf_qtr, hf_qrs, hf_qst; bool qtr_ok = false; bool qrs_ok = false; bool qst_ok = false;
	for(unsigned i=0;i<hfh_vec.size();++i)
	{
		if( hfh_vec[i] != hfh1 )
		{
			OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfh_vec[i]); int v_id0 = hfv_it->idx();
			++hfv_it; int v_id1 = hfv_it->idx(); ++hfv_it; int v_id2 = hfv_it->idx();
			if( !qtr_ok )
			{
				if( (v_id0 == q_id && v_id1 == t_id && v_id2 == r_id)
					|| (v_id1 == q_id && v_id2 == t_id && v_id0 == r_id)
					|| (v_id2 == q_id && v_id0 == t_id && v_id1 == r_id) )
				{
					hf_qtr = hfh_vec[i];
					qtr_ok = true;
				}
			}

			if( !qrs_ok )
			{
				if( (v_id0 == q_id && v_id1 == r_id && v_id2 == s_id)
					|| (v_id1 == q_id && v_id2 == r_id && v_id0 == s_id)
					|| (v_id2 == q_id && v_id0 == r_id && v_id1 == s_id) )
				{
					hf_qrs = hfh_vec[i];
					qrs_ok = true;
				}
			}

			if( !qst_ok )
			{
				if( (v_id0 == q_id && v_id1 == s_id && v_id2 == t_id)
					|| (v_id1 == q_id && v_id2 == s_id && v_id0 == t_id)
					|| (v_id2 == q_id && v_id0 == s_id && v_id1 == t_id) )
				{
					hf_qst = hfh_vec[i];
					qst_ok = true;
				}
			}
		}
	}

	if( ch1.idx() > ch0.idx() )
	{
		mesh_->delete_cell(ch1); mesh_->delete_cell(ch0);
	}
	else
	{
		mesh_->delete_cell(ch0); mesh_->delete_cell(ch1);
	}
	

	std::vector<OpenVolumeMesh::VertexHandle > one_face(3);
	one_face[0] = p_vh; one_face[1] = q_vh;
	one_face[2] = OpenVolumeMesh::VertexHandle(r_id);
	OpenVolumeMesh::FaceHandle f_pqr = mesh_->add_face(one_face);
	OpenVolumeMesh::HalfFaceHandle hf_pqr = mesh_->halfface_handle(f_pqr, 0);

	one_face[2] = OpenVolumeMesh::VertexHandle(s_id);
	OpenVolumeMesh::FaceHandle f_pqs = mesh_->add_face(one_face);
	OpenVolumeMesh::HalfFaceHandle hf_pqs = mesh_->halfface_handle(f_pqs, 0);

	one_face[2] = OpenVolumeMesh::VertexHandle(t_id);
	OpenVolumeMesh::FaceHandle f_pqt = mesh_->add_face(one_face);
	OpenVolumeMesh::HalfFaceHandle hf_pqt = mesh_->halfface_handle(f_pqt, 0);

	std::vector<OpenVolumeMesh::HalfFaceHandle> one_cell(4);
	one_cell[0] = hf_prt; one_cell[1] = hf_qtr; 
	one_cell[2] = hf_pqr; one_cell[3] = mesh_->opposite_halfface_handle(hf_pqt);
	mesh_->add_cell(one_cell);

	one_cell[0] = hf_psr; one_cell[1] = hf_qrs; 
	one_cell[2] = hf_pqs; one_cell[3] = mesh_->opposite_halfface_handle(hf_pqr);
	mesh_->add_cell(one_cell);

	one_cell[0] = hf_pts; one_cell[1] = hf_qst; 
	one_cell[2] = hf_pqt; one_cell[3] = mesh_->opposite_halfface_handle(hf_pqs);
	mesh_->add_cell(one_cell);

	mesh_->delete_face(fh);

	return true;
}

bool is_flip_edge_32_ok(VolumeMesh* mesh_, OpenVolumeMesh::EdgeHandle eh)
{
	if(mesh_->is_boundary(eh) || mesh_->valence(eh) != 3)
		return false;

	OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge(eh);
	OpenVolumeMesh::VertexHandle p_vh = edge.to_vertex();   OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(p_vh);
	OpenVolumeMesh::VertexHandle q_vh = edge.from_vertex(); OpenVolumeMesh::Geometry::Vec3d q = mesh_->vertex(q_vh);

	OpenVolumeMesh::HalfEdgeHandle heh = mesh_->halfedge_handle(eh, 0);
	OpenVolumeMesh::HalfEdgeHalfFaceIter hehf_it = mesh_->hehf_iter(heh); 
	OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hehf_it.cur_handle());
	OpenVolumeMesh::VertexHandle tri_v[3];
	tri_v[0] = hfv_it.cur_handle();
	++hfv_it; tri_v[1] = hfv_it.cur_handle();
	++hfv_it; tri_v[2] = hfv_it.cur_handle();
	
	if( (tri_v[0] == q_vh && tri_v[1] == p_vh) 
		|| (tri_v[1] == q_vh && tri_v[2] == p_vh)
		|| (tri_v[2] == q_vh && tri_v[0] == p_vh) )
	{}
	else
	{
		heh = mesh_->opposite_halfedge_handle(heh);
	}
	
	int i = 0; OpenVolumeMesh::Geometry::Vec3d tri_p[3];
	for( hehf_it = mesh_->hehf_iter(heh);hehf_it; ++hehf_it )
	{
		for( hfv_it = mesh_->hfv_iter(hehf_it.cur_handle()); hfv_it; ++hfv_it )
		{
			if( hfv_it.cur_handle() != p_vh && hfv_it.cur_handle() != q_vh)
			{
				tri_p[i] = mesh_->vertex( hfv_it.cur_handle() ); ++i; 
				break;
			}
		}
	}

	if( !line_intersect_triangle(p, q, tri_p[0], tri_p[1], tri_p[2]) )
	{
		return false;
	}
	
	return true;
}

bool flip_edge_32(VolumeMesh* mesh_, OpenVolumeMesh::EdgeHandle eh)
{
	if( !is_flip_edge_32_ok(mesh_, eh) ) 
		return false;

	OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge(eh);
	OpenVolumeMesh::VertexHandle p_vh = edge.to_vertex();   OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(p_vh);
	OpenVolumeMesh::VertexHandle q_vh = edge.from_vertex(); OpenVolumeMesh::Geometry::Vec3d q = mesh_->vertex(q_vh);

	OpenVolumeMesh::HalfEdgeHandle heh = mesh_->halfedge_handle(eh, 0);
	OpenVolumeMesh::HalfEdgeHalfFaceIter hehf_it = mesh_->hehf_iter(heh); 
	OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hehf_it.cur_handle());
	std::vector<OpenVolumeMesh::VertexHandle> tri_v(3);
	tri_v[0] = hfv_it.cur_handle();
	++hfv_it; tri_v[1] = hfv_it.cur_handle();
	++hfv_it; tri_v[2] = hfv_it.cur_handle();

	if( (tri_v[0] == q_vh && tri_v[1] == p_vh) 
		|| (tri_v[1] == q_vh && tri_v[2] == p_vh)
		|| (tri_v[2] == q_vh && tri_v[0] == p_vh) )
	{}
	else
	{
		heh = mesh_->opposite_halfedge_handle(heh);
	}

	int i = 0; OpenVolumeMesh::CellHandle ch[3]; OpenVolumeMesh::FaceHandle fh[3];
	std::vector<OpenVolumeMesh::HalfFaceHandle> p_one_cell; p_one_cell.reserve(4);
	std::vector<OpenVolumeMesh::HalfFaceHandle> q_one_cell; q_one_cell.reserve(4);
	for( hehf_it = mesh_->hehf_iter(heh);hehf_it; ++hehf_it )
	{
		ch[i] = mesh_->incident_cell( hehf_it.cur_handle() );
		fh[i] = mesh_->face_handle(hehf_it.cur_handle() );
		OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(ch[i]);
		std::vector<OpenVolumeMesh::HalfFaceHandle> hfh_vec = cell.halffaces();
		for( unsigned j=0; j<hfh_vec.size(); ++j)
		{
			bool p_ok = false; bool q_ok = false;
			for( hfv_it = mesh_->hfv_iter(hfh_vec[j]); hfv_it; ++hfv_it)
			{
				if(hfv_it.cur_handle() == p_vh )
				{
					p_ok = true;
				}
				else if(hfv_it.cur_handle() == q_vh)
				{
					q_ok = true;
				}
			}

			if(p_ok && !q_ok)
			{
				p_one_cell.push_back(hfh_vec[j]);
			}
			else if(!p_ok && q_ok)
			{
				q_one_cell.push_back(hfh_vec[j]);
			}
		}

		for( hfv_it = mesh_->hfv_iter(hehf_it.cur_handle()); hfv_it; ++hfv_it )
		{
			if( hfv_it.cur_handle() != p_vh && hfv_it.cur_handle() != q_vh)
			{
				tri_v[i] = hfv_it.cur_handle(); ++i; 
				break;
			}
		}
	}

	int delete_index[3];
	if( ch[0].idx() > ch[1].idx() )
	{
		if(ch[0].idx() > ch[2].idx())
		{
			delete_index[0] = 0;
			if( ch[1].idx() > ch[2].idx() )
			{
				delete_index[1] = 1; delete_index[2] = 2;
			}
			else
			{
				delete_index[1] = 2; delete_index[2] = 1;
			}
		}
		else
		{
			delete_index[0] = 2; delete_index[1] = 0; delete_index[2] = 1;
		}
	}
	else
	{
		if(ch[1].idx() > ch[2].idx())
		{
			delete_index[0] = 1;
			if( ch[0].idx() > ch[2].idx() )
			{
				delete_index[1] = 0; delete_index[2] = 2;
			}
			else
			{
				delete_index[1] = 2; delete_index[2] = 0;
			}
		}
		else
		{
			delete_index[0] = 2; delete_index[1] = 1; delete_index[2] = 0;
		}
	}

	mesh_->delete_cell( ch[delete_index[0]] ); mesh_->delete_cell( ch[delete_index[1]] ); mesh_->delete_cell( ch[delete_index[2]] );

	OpenVolumeMesh::FaceHandle tfh = mesh_->add_face( tri_v );
	p_one_cell.push_back( mesh_->halfface_handle(tfh, 0) );
	q_one_cell.push_back( mesh_->halfface_handle(tfh, 1) );

	mesh_->add_cell(p_one_cell); mesh_->add_cell(q_one_cell);

	if( fh[0].idx() > fh[1].idx() )
	{
		if(fh[0].idx() > fh[2].idx())
		{
			delete_index[0] = 0;
			if( fh[1].idx() > fh[2].idx() )
			{
				delete_index[1] = 1; delete_index[2] = 2;
			}
			else
			{
				delete_index[1] = 2; delete_index[2] = 1;
			}
		}
		else
		{
			delete_index[0] = 2; delete_index[1] = 0; delete_index[2] = 1;
		}
	}
	else
	{
		if(fh[1].idx() > fh[2].idx())
		{
			delete_index[0] = 1;
			if( fh[0].idx() > fh[2].idx() )
			{
				delete_index[1] = 0; delete_index[2] = 2;
			}
			else
			{
				delete_index[1] = 2; delete_index[2] = 0;
			}
		}
		else
		{
			delete_index[0] = 2; delete_index[1] = 1; delete_index[2] = 0;
		}
	}
	mesh_->delete_face(fh[delete_index[0]]); mesh_->delete_face(fh[delete_index[1]]); mesh_->delete_face(fh[delete_index[2]]);

	mesh_->delete_edge(eh);
	return true;
}

bool baryCoord( const OpenMesh::Vec3d& _p, const OpenMesh::Vec3d& _u, const OpenMesh::Vec3d& _v, const OpenMesh::Vec3d& _w, OpenMesh::Vec3d&_result )
{
	OpenMesh::Vec3d  vu = _v - _u, wu = _w - _u, pu = _p - _u;

	// find largest absolute coodinate of normal
	double nx = vu[1]*wu[2] - vu[2]*wu[1],
		ny = vu[2]*wu[0] - vu[0]*wu[2],
		nz = vu[0]*wu[1] - vu[1]*wu[0],
		ax = fabs(nx),
		ay = fabs(ny),
		az = fabs(nz);


	unsigned char max_coord;

	if ( ax > ay )
	{
		if ( ax > az ) 
		{
			max_coord = 0;
		}
		else
		{
			max_coord = 2;
		}
	}
	else
	{
		if ( ay > az )
		{
			max_coord = 1;
		}
		else 
		{
			max_coord = 2;
		}
	}


	// solve 2D problem
	switch (max_coord)
	{
	case 0:
		{      
			if (1.0+ax == 1.0) return false;
			_result[1] = 1.0 + (pu[1]*wu[2]-pu[2]*wu[1])/nx - 1.0;
			_result[2] = 1.0 + (vu[1]*pu[2]-vu[2]*pu[1])/nx - 1.0;
			_result[0] = 1.0 - _result[1] - _result[2];
		}
		break;

	case 1:
		{
			if (1.0+ay == 1.0) return false;
			_result[1] = 1.0 + (pu[2]*wu[0]-pu[0]*wu[2])/ny - 1.0;
			_result[2] = 1.0 + (vu[2]*pu[0]-vu[0]*pu[2])/ny - 1.0;
			_result[0] = 1.0 - _result[1] - _result[2];
		}
		break;

	case 2:
		{
			if (1.0+az == 1.0) return false;
			_result[1] = 1.0 + (pu[0]*wu[1]-pu[1]*wu[0])/nz - 1.0;
			_result[2] = 1.0 + (vu[0]*pu[1]-vu[1]*pu[0])/nz - 1.0;
			_result[0] = 1.0 - _result[1] - _result[2];
		}
		break;
	}

	return true;
}

double distPointTriangleSquared(const OpenVolumeMesh::Geometry::Vec3d& _p,
	const OpenVolumeMesh::Geometry::Vec3d& _v0,
	const OpenVolumeMesh::Geometry::Vec3d& _v1,
	const OpenVolumeMesh::Geometry::Vec3d& _v2,
	OpenVolumeMesh::Geometry::Vec3d& _nearestPoint)
{
	OpenVolumeMesh::Geometry::Vec3d v0v1(_v1[0] - _v0[0], _v1[1] - _v0[1], _v1[2] - _v0[2]);
	OpenVolumeMesh::Geometry::Vec3d v0v2 = _v2 - _v0;
	OpenVolumeMesh::Geometry::Vec3d n = v0v1 % v0v2; // not normalized !
	double d = n.sqrnorm();

	// Check if the triangle is degenerated
	if (d < FLT_MIN && d > -FLT_MIN)
	{
		std::cerr << "distPointTriangleSquared: Degenerated triangle !\n";
		return -1.0;
	}
	double invD = 1.0 / d;

	// these are not needed for every point, should still perform
	// better with many points against one triangle
	OpenVolumeMesh::Geometry::Vec3d v1v2 = _v2 - _v1;
	double inv_v0v2_2 = 1.0 / v0v2.sqrnorm();
	double inv_v0v1_2 = 1.0 / v0v1.sqrnorm();
	double inv_v1v2_2 = 1.0 / v1v2.sqrnorm();

	OpenVolumeMesh::Geometry::Vec3d v0p = _p - _v0;
	OpenVolumeMesh::Geometry::Vec3d t = v0p % n;
	double  s01, s02, s12;
	double a = (t | v0v2) * -invD;
	double b = (t | v0v1) * invD;

	if (a < 0)
	{
		// Calculate the distance to an edge or a corner vertex
		s02 = (v0v2 | v0p) * inv_v0v2_2;
		if (s02 < 0.0)
		{
			s01 = (v0v1 | v0p) * inv_v0v1_2;
			if (s01 <= 0.0)
			{
				v0p = _v0;
			}
			else if (s01 >= 1.0)
			{
				v0p = _v1;
			}
			else
			{
				v0p = _v0 + v0v1 * s01;
			}
		}
		else if (s02 > 1.0)
		{
			s12 = (v1v2 | (_p - _v1)) * inv_v1v2_2;
			if (s12 >= 1.0)
			{
				v0p = _v2;
			}
			else if (s12 <= 0.0)
			{
				v0p = _v1;
			}
			else
			{
				v0p = _v1 + v1v2 * s12;
			}
		}
		else
		{
			v0p = _v0 + v0v2 * s02;
		}
	}
	else if (b < 0.0)
	{
		// Calculate the distance to an edge or a corner vertex
		s01 = (v0v1 | v0p) * inv_v0v1_2;
		if (s01 < 0.0)
		{
			s02 = (v0v2 | v0p) * inv_v0v2_2;
			if (s02 <= 0.0)
			{
				v0p = _v0;
			}
			else if (s02 >= 1.0)
			{
				v0p = _v2;
			}
			else
			{
				v0p = _v0 + v0v2 * s02;
			}
		}
		else if (s01 > 1.0)
		{
			s12 = (v1v2 | (_p - _v1)) * inv_v1v2_2;
			if (s12 >= 1.0)
			{
				v0p = _v2;
			}
			else if (s12 <= 0.0)
			{
				v0p = _v1;
			}
			else
			{
				v0p = _v1 + v1v2 * s12;
			}
		}
		else
		{
			v0p = _v0 + v0v1 * s01;
		}
	}
	else if (a + b > 1.0)
	{
		// Calculate the distance to an edge or a corner vertex
		s12 = (v1v2 | (_p - _v1)) * inv_v1v2_2;
		if (s12 >= 1.0)
		{
			s02 = (v0v2 | v0p) * inv_v0v2_2;
			if (s02 <= 0.0)
			{
				v0p = _v0;
			}
			else if (s02 >= 1.0)
			{
				v0p = _v2;
			}
			else
			{
				v0p = _v0 + v0v2*s02;
			}
		}
		else if (s12 <= 0.0)
		{
			s01 = (v0v1 | v0p) * inv_v0v1_2;
			if (s01 <= 0.0)
			{
				v0p = _v0;
			}
			else if (s01 >= 1.0)
			{
				v0p = _v1;
			}
			else
			{
				v0p = _v0 + v0v1 * s01;
			}
		}
		else
		{
			v0p = _v1 + v1v2 * s12;
		}
	}
	else
	{
		// Calculate the distance to an interior point of the triangle
		_nearestPoint = _p - n*((n | v0p) * invD);
		return (_nearestPoint - _p).sqrnorm();
	}

	_nearestPoint = v0p;

	return (_nearestPoint - _p).sqrnorm();
}

double distPointTriangleSquared2(const double _p[3], const double _v0[3], const double _v1[3], const double _v2[3], double _nearestPoint[3])
{
	double v0v1[3] = { _v1[0] - _v0[0], _v1[1] - _v0[1], _v1[2] - _v0[2] };
	double v0v2[3] = { _v2[0] - _v0[0], _v2[1] - _v0[1], _v2[2] - _v0[2] };
	double n[3] = { v0v1[1] * v0v2[2] - v0v1[2] * v0v2[1], v0v1[2] * v0v2[0] - v0v1[0] * v0v2[2], v0v1[0] * v0v2[1] - v0v1[1] * v0v2[0] };
	double d = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
	/*OpenVolumeMesh::Geometry::Vec3d v0v1(_v1[0] - _v0[0], _v1[1] - _v0[1], _v1[2] - _v0[2]);
	OpenVolumeMesh::Geometry::Vec3d v0v2 = _v2 - _v0;
	OpenVolumeMesh::Geometry::Vec3d n = v0v1 % v0v2; // not normalized !
	double d = n.sqrnorm();*/

	// Check if the triangle is degenerated
	if (d < FLT_MIN && d > -FLT_MIN)
	{
		std::cerr << "distPointTriangleSquared: Degenerated triangle !\n";
		return -1.0;
	}
	double invD = 1.0 / d;

	// these are not needed for every point, should still perform
	// better with many points against one triangle
	double v1v2[3] = { _v2[0] - _v1[0], _v2[1] - _v1[1], _v2[2] - _v1[2] };
	double inv_v0v2_2 = 1.0 / (v0v2[0] * v0v2[0] + v0v2[1] * v0v2[1] + v0v2[2] * v0v2[2]);
	double inv_v0v1_2 = 1.0 / (v0v1[0] * v0v1[0] + v0v1[1] * v0v1[1] + v0v1[2] * v0v1[2]);
	double inv_v1v2_2 = 1.0 / (v1v2[0] * v1v2[0] + v1v2[1] * v1v2[1] + v1v2[2] * v1v2[2]);
	double v0p[3] = { _p[0] - _v0[0], _p[1] - _v0[1], _p[2] - _v0[2] };
	double t[3] = { v0p[1] * n[2] - v0p[2] * n[1], v0p[2] * n[0] - v0p[0] * n[2], v0p[0] * n[1] - v0p[1] * n[0] };
	double  s01, s02, s12;
	double a = (t[0] * v0v2[0] + t[1] * v0v2[1] + t[2] * v0v2[2]) * -invD;
	double b = (t[0] * v0v1[0] + t[1] * v0v1[1] + t[2] * v0v1[2]) * invD;

	/*OpenVolumeMesh::Geometry::Vec3d v1v2 = _v2 - _v1;
	double inv_v0v2_2 = 1.0 / v0v2.sqrnorm();
	double inv_v0v1_2 = 1.0 / v0v1.sqrnorm();
	double inv_v1v2_2 = 1.0 / v1v2.sqrnorm();

	OpenVolumeMesh::Geometry::Vec3d v0p = _p - _v0;
	OpenVolumeMesh::Geometry::Vec3d t = v0p % n;
	double  s01, s02, s12;
	double a = (t | v0v2) * -invD;
	double b = (t | v0v1) * invD;*/

	if (a < 0)
	{
		// Calculate the distance to an edge or a corner vertex
		s02 = (v0p[0] * v0v2[0] + v0p[1] * v0v2[1] + v0p[2] * v0v2[2]) * inv_v0v2_2;
		if (s02 < 0.0)
		{
			s01 = (v0p[0] * v0v1[0] + v0p[1] * v0v1[1] + v0p[2] * v0v1[2]) * inv_v0v1_2;
			if (s01 <= 0.0)
			{
				//v0p = _v0;
				v0p[0] = _v0[0]; v0p[1] = _v0[1]; v0p[2] = _v0[2];
			}
			else if (s01 >= 1.0)
			{
				//v0p = _v1;
				v0p[0] = _v1[0]; v0p[1] = _v1[1]; v0p[2] = _v1[2];
			}
			else
			{
				//v0p = _v0 + v0v1 * s01;
				v0p[0] = _v0[0] + v0v1[0] * s01; v0p[1] = _v0[1] + v0v1[1] * s01; v0p[2] = _v0[2] + v0v1[2] * s01;
			}
		}
		else if (s02 > 1.0)
		{
			s12 = (v1v2[0] * (_p[0] - _v1[0]) + v1v2[1] * (_p[1] - _v1[1]) + v1v2[2] * (_p[2] - _v1[2])) * inv_v1v2_2;
			if (s12 >= 1.0)
			{
				//v0p = _v2;
				v0p[0] = _v2[0]; v0p[1] = _v2[1]; v0p[2] = _v2[2];
			}
			else if (s12 <= 0.0)
			{
				//v0p = _v1;
				v0p[0] = _v1[0]; v0p[1] = _v1[1]; v0p[2] = _v1[2];
			}
			else
			{
				//v0p = _v1 + v1v2 * s12;
				v0p[0] = _v1[0] + v1v2[0] * s12; v0p[1] = _v1[1] + v1v2[1] * s12; v0p[2] = _v1[2] + v1v2[2] * s12;
			}
		}
		else
		{
			//v0p = _v0 + v0v2 * s02;
			v0p[0] = _v0[0] + v0v2[0] * s02; v0p[1] = _v0[1] + v0v2[1] * s02; v0p[2] = _v0[2] + v0v2[2] * s02;
		}
	}
	else if (b < 0.0)
	{
		// Calculate the distance to an edge or a corner vertex
		//s01 = (v0v1 | v0p) * inv_v0v1_2;
		s01 = (v0v1[0] * v0p[0] + v0v1[1] * v0p[1] + v0v1[2] * v0p[2]) * inv_v0v1_2;
		if (s01 < 0.0)
		{
			//s02 = (v0v2 | v0p) * inv_v0v2_2;
			s02 = (v0v2[0] * v0p[0] + v0v2[1] * v0p[1] + v0v2[2] * v0p[2]) * inv_v0v2_2;
			if (s02 <= 0.0)
			{
				//v0p = _v0;
				v0p[0] = _v0[0]; v0p[1] = _v0[1]; v0p[2] = _v0[2];
			}
			else if (s02 >= 1.0)
			{
				//v0p = _v2;
				v0p[0] = _v2[0]; v0p[1] = _v2[1]; v0p[2] = _v2[2];
			}
			else
			{
				//v0p = _v0 + v0v2 * s02;
				v0p[0] = _v0[0] + v0v2[0] * s02; v0p[1] = _v0[1] + v0v2[1] * s02; v0p[2] = _v0[2] + v0v2[2] * s02;
			}
		}
		else if (s01 > 1.0)
		{
			//s12 = (v1v2 | (_p - _v1)) * inv_v1v2_2;
			s12 = (v1v2[0] * (_p[0] - _v1[0]) + v1v2[1] * (_p[1] - _v1[1]) + v1v2[2] * (_p[2] - _v1[2])) * inv_v1v2_2;
			if (s12 >= 1.0)
			{
				//v0p = _v2;
				v0p[0] = _v2[0]; v0p[1] = _v2[1]; v0p[2] = _v2[2];
			}
			else if (s12 <= 0.0)
			{
				//v0p = _v1;
				v0p[0] = _v1[0]; v0p[1] = _v1[1]; v0p[2] = _v1[2];
			}
			else
			{
				//v0p = _v1 + v1v2 * s12;
				v0p[0] = _v1[0] + v1v2[0] * s12; v0p[1] = _v1[1] + v1v2[1] * s12; v0p[2] = _v1[2] + v1v2[2] * s12;
			}
		}
		else
		{
			//v0p = _v0 + v0v1 * s01;
			v0p[0] = _v0[0] + v0v1[0] * s01; v0p[1] = _v0[1] + v0v1[1] * s01; v0p[2] = _v0[2] + v0v1[2] * s01;
		}
	}
	else if (a + b > 1.0)
	{
		// Calculate the distance to an edge or a corner vertex
		//s12 = (v1v2 | (_p - _v1)) * inv_v1v2_2;
		s12 = (v1v2[0] * (_p[0] - _v1[0]) + v1v2[1] * (_p[1] - _v1[1]) + v1v2[2] * (_p[2] - _v1[2])) * inv_v1v2_2;
		if (s12 >= 1.0)
		{
			//s02 = (v0v2 | v0p) * inv_v0v2_2;
			s02 = (v0v2[0] * v0p[0] + v0v2[1] * v0p[1] + v0v2[2] * v0p[2]) * inv_v0v2_2;
			if (s02 <= 0.0)
			{
				//v0p = _v0;
				v0p[0] = _v0[0]; v0p[1] = _v0[1]; v0p[2] = _v0[2];
			}
			else if (s02 >= 1.0)
			{
				//v0p = _v2;
				v0p[0] = _v2[0]; v0p[1] = _v2[1]; v0p[2] = _v2[2];
			}
			else
			{
				//v0p = _v0 + v0v2*s02;
				v0p[0] = _v0[0] + v0v2[0] * s02; v0p[1] = _v0[1] + v0v2[1] * s02; v0p[2] = _v0[2] + v0v2[2] * s02;
			}
		}
		else if (s12 <= 0.0)
		{
			//s01 = (v0v1 | v0p) * inv_v0v1_2;
			s01 = (v0v1[0] * v0p[0] + v0v1[1] * v0p[1] + v0v1[2] * v0p[2]) * inv_v0v1_2;
			if (s01 <= 0.0)
			{
				//v0p = _v0;
				v0p[0] = _v0[0]; v0p[1] = _v0[1]; v0p[2] = _v0[2];
			}
			else if (s01 >= 1.0)
			{
				//v0p = _v1;
				v0p[0] = _v1[0]; v0p[1] = _v1[1]; v0p[2] = _v1[2];
			}
			else
			{
				//v0p = _v0 + v0v1 * s01;
				v0p[0] = _v0[0] + v0v1[0] * s01; v0p[1] = _v0[1] + v0v1[1] * s01; v0p[2] = _v0[2] + v0v1[2] * s01;
			}
		}
		else
		{
			//v0p = _v1 + v1v2 * s12;
			v0p[0] = _v1[0] + v1v2[0] * s12; v0p[1] = _v1[1] + v1v2[1] * s12; v0p[2] = _v1[2] + v1v2[2] * s12;
		}
	}
	else
	{
		// Calculate the distance to an interior point of the triangle
		double temp = (n[0] * v0p[0] + n[1] * v0p[1] + n[2] * v0p[2])*invD;
		_nearestPoint[0] = _p[0] - n[0] * temp; _nearestPoint[1] = _p[1] - n[1] * temp; _nearestPoint[2] = _p[2] - n[2] * temp;
		//return (_nearestPoint - _p).sqrnorm();
		return ((_nearestPoint[0] - _p[0])*(_nearestPoint[0] - _p[0]) + (_nearestPoint[1] - _p[1])*(_nearestPoint[1] - _p[1]) + (_nearestPoint[2] - _p[2])*(_nearestPoint[2] - _p[2]));
	}

	_nearestPoint[0] = v0p[0]; _nearestPoint[1] = v0p[1]; _nearestPoint[2] = v0p[2];

	return ((_nearestPoint[0] - _p[0])*(_nearestPoint[0] - _p[0]) + (_nearestPoint[1] - _p[1])*(_nearestPoint[1] - _p[1]) + (_nearestPoint[2] - _p[2])*(_nearestPoint[2] - _p[2]));
}