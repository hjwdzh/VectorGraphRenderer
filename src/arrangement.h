#ifndef VECTORGRAPH_RENDER_ARRANGEMENT_H_
#define VECTORGRAPH_RENDER_ARRANGEMENT_H_

#include <vector>

#include <Eigen/Core>

#include "types.h"

void ComputeTriangleArrangement(const std::vector<Eigen::Vector3d>& vertices, const Eigen::Vector3i& f, int fid, Arrangement_2* pout);

void MergeFaceInsideArrangement(Arrangement_2& out);

#endif