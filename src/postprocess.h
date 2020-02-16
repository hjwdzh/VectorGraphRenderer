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
		std::vector<int> outer_type;
		std::vector<std::vector<int> > inner_indices;
		std::vector<std::vector<int> > inner_type;
	};

	void CollectFaceAndVertices(const Mesh& mesh, const Arrangement_2& overlay, double angle_thres);
	void RemoveRedundantVertices();
	void MergeDuplex(const Mesh& mesh);
	void CollectEdges(const Mesh& mesh);

	void SaveToFile(const Mesh& mesh, const char* filename);
	void SaveToSVG(const Mesh& mesh, const char* filename);
protected:
	void RemoveRedundantVerticesFromLoop(std::vector<int>& indices, std::vector<int>& types);
	void ComputeDegree();

private:
	std::vector<EdgeList> facets_;

	std::vector<int> facets_id_;
	std::vector<VertexSignature > points_;
	std::vector<K> depths_;

	std::vector<int> degrees_;

	std::map<std::pair<int,int>, std::set<int> > edge2face_;
};

#endif