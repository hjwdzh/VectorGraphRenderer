#ifndef VECTORGRAPH_RENDER_ARRANGEMENT_H_
#define VECTORGRAPH_RENDER_ARRANGEMENT_H_

#include <vector>

#include <Eigen/Core>

#include "mesh.h"
#include "types.h"

void ComputeTriangleArrangement(const Mesh& mesh, int fid, Arrangement_2* pout);

void MergeFaceInsideArrangement(Arrangement_2& out);

#endif