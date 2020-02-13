#ifndef VECTORGRAPH_RENDER_POSTPROCESS_H_
#define VECTORGRAPH_RENDER_POSTPROCESS_H_

#include "mesh.h"
#include "types.h"

class PostProcess
{
public:
	typedef std::pair<Arrangement_2::Vertex_const_handle,int> VertexSignature;
	struct EdgeList {
		std::vector<int> outer_indices;
		std::vector<std::vector<int> > inner_indices;
	};

	void CollectFaceAndVertices(const Arrangement_2& overlay);
	void RemoveRedundantVertices();
	void MergeDuplex(const Mesh& mesh);

	void SaveToFile(const Mesh& mesh, double angle_thres, const char* filename);
protected:
	void RemoveRedundantVerticesFromLoop(std::vector<int>& indices);
	void ComputeDegree();

private:
	std::vector<EdgeList> facets_;

	std::vector<int> facets_id_;
	std::vector<VertexSignature > points_;

	std::vector<int> degrees_;
};

#endif