#ifndef VECTORGRAPH_RENDERER_CAMERA_H_
#define VECTORGRAPH_RENDERER_CAMERA_H_

#include <Eigen/Core>
class Mesh;
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
	void Undistort(Eigen::Vector3d& v) const {
		v[0] = (v[0] * width_ - cx_) / fx_ * v[2];
		v[1] = (v[1] * width_ - cx_) / fx_ * v[2];
	}
private:
	float fx_, fy_, cx_, cy_;
	int width_, height_;
	Eigen::Matrix4d world2cam_;

	double angle_;
};

#endif
