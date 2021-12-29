#pragma once
#include <Eigen/core>
#include <Eigen/Geometry>
#include <Eigen/dense>


class Movable
{
public:
	Movable();
	Movable(const Movable& mov);
	Eigen::Matrix4f MakeTransScale();
	Eigen::Matrix4d MakeTransd();
	Eigen::Matrix4d MakeTransScaled();
	void MyTranslate(Eigen::Vector3d amt, bool preRotation);



	// ASSIGNMENT 1 TASK 1 - Rotate annoucnment
	Eigen::Matrix3d GetRotation() { return Tout.rotation(); }

	void TranslateInSystem(Eigen::Matrix3d , Eigen::Vector3d amt);

	 
	void RotateInSystem(Eigen::Vector3d rotAxis, double angle);
	void MyRotate(Eigen::Vector3f rotAxis, float angle);
	// !!!!!!!!!!!!

	void MyRotate(Eigen::Vector3d rotAxis, double angle);
	void MyRotate(const Eigen::Matrix3d &rot);
	void MyScale(Eigen::Vector3d amt);

	Eigen::Matrix3d GetRotation() const{ return Tout.rotation().matrix(); }

	virtual ~Movable() {}

	Eigen::Affine3d Tout,Tin;

private:
};

