#pragma once

#include <stdio.h>
#include <vector>
#include <algorithm>
#include "Mesh.h"
#include "AABB.h"

using namespace std;

class SurfaceRemeshing
{
public:
	SurfaceRemeshing(void);
	~SurfaceRemeshing(void);
	SurfaceRemeshing(const char *subject, const char *sphere, const char *dfield, bool keepColor, const char *sphere_t = NULL, const char *colormap = NULL, vector<string> property = vector<string>(), bool interpolation = true, bool backward = false);
	void saveDeformedSurface(const char *filename);
	void saveDeformedSphere(const char *filename);
	void saveDeformedProperty(const char *filename, bool header = true);
	
private:
	void reconsCoord(const float *v0, float *v1, float *Y, float *coeff, int degree, float *pole);
	float dataInterpolation(float *refMap, int index, float *coeff, Mesh *mesh);
	int dataInterpolation(vector<int *> refMap, int index, float *coeff, Mesh *mesh, int channel);
	float dataMedian(float *refMap, int index, Mesh *mesh);
	void deformSurface(void);
	void deformData(void);

private:
	int m_degree;
	float m_pole[3];
	float *m_coeff;
	float *m_x, *m_y, *m_z;
	AABB *m_tree;
	bool m_interpolation;
	bool m_backward;
	bool m_keepColor;
	Mesh *m_sphere, *m_sphere_subj, *m_subj, *m_remesh;
	vector<float *> m_refMap;
	vector<float *> m_deData;
	vector<string> m_property;
	vector<int *> m_color;
	vector<int *> m_color_base;
};


