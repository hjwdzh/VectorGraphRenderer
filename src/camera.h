#ifndef VECTORGRAPH_RENDERER_CAMERA_H_
#define VECTORGRAPH_RENDERER_CAMERA_H_

#include "mesh.h"
#include "types.h"

class Camera
{
public:
	Camera();
	void LoadFromFile(const char* filename);
	void ApplyExtrinsic(Mesh& mesh);
	void ApplyIntrinsic(Mesh& mesh);

	double GetAngle() const {
		return angle_;
	}
private:
	float fx_, fy_, cx_, cy_;
	int width_, height_;
	Eigen::Matrix4d world2cam_;

	double angle_;
};

#endif
