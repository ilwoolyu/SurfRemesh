#pragma once
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <float.h>
#include "Mesh.h"
#include "Geom.h"
#include "string.h"

using namespace std;

class AABB
{
public:
	AABB(void);
	AABB(Mesh *mesh)
	{
		m_mesh = mesh;
		initTree();
		update();
	}
	~AABB(void);
	int closestFace(float *v, float *coeff, float range = 0)
	{
		float err = 0; // numerical error
		bool success = false;
		int index = 0;

		float vmin[3] = {v[0] - range, v[1] - range, v[2] - range};
		float vmax[3] = {v[0] + range, v[1] + range, v[2] + range};
		vector<int> cand;
		float eps = 0;
		for (int trial = 0; cand.empty(); trial++)
		{
			searchTree(vmin, vmax, m_tree, &cand, eps);
			if (trial == 0) eps = 1e-7;
			else eps *= 10;
		}
		sort(cand.begin(), cand.end());
		cand.erase(unique(cand.begin(), cand.end()), cand.end());

		for (int trial = 0; !success; trial++)
		{
			float min_dist = FLT_MAX;
			for (int i = 0; i < cand.size(); i++)
			{
				const float *a = m_mesh->face(cand[i])->vertex(0)->fv();
				const float *b = m_mesh->face(cand[i])->vertex(1)->fv();
				const float *c = m_mesh->face(cand[i])->vertex(2)->fv();

				// closest distance
				float dist = Coordinate::dpoint2tri(a, b, c, v);

				if (dist < min_dist)
				{
					index = cand[i];
					min_dist = dist;
				}
			}
			success = true;
		}

		Face f = *m_mesh->face(index);
		const Vertex &a = *f.vertex(0);
		const Vertex &b = *f.vertex(1);
		const Vertex &c = *f.vertex(2);

		// bary centric
		Coordinate::cart2bary((float *)a.fv(), (float *)b.fv(), (float *)c.fv(), v, coeff);

		return index;
	}
	void update()
	{
		updateTree(m_tree);
	}

private:
	struct node
	{
		node *left;
		node *right;
		float x0, y0, z0;
		float x1, y1, z1;
		vector<int> cand;
	};
	void initTree(void)
	{
		vector<int> cand;
		vector<float *> range;
		for (int i = 0; i < m_mesh->nFace(); i++)
		{
			Face f = *m_mesh->face(i);
			float *r = new float[6];
			boundingBox(f, r);
			range.push_back(r);
		}
		for (int i = 0; i < m_mesh->nFace(); i++) cand.push_back(i);

		m_tree = construction(range, cand);
	};
	void searchTree(float *pmin, float *pmax, node *root, vector<int> *cand, float eps = 0, bool trace = false)
	{
		bool traceOverL = false, traceOverR = false;

		if (!trace)
		{
			if (root->left != NULL && 
				pmin[0] >= root->left->x0 - eps && pmin[1] >= root->left->y0 - eps && pmin[2] >= root->left->z0 - eps && 
				pmin[0] <= root->left->x1 + eps && pmin[1] <= root->left->y1 + eps && pmin[2] <= root->left->z1 + eps &&
				pmax[0] >= root->left->x0 - eps && pmax[1] >= root->left->y0 - eps && pmax[2] >= root->left->z0 - eps && 
				pmax[0] <= root->left->x1 + eps && pmax[1] <= root->left->y1 + eps && pmax[2] <= root->left->z1 + eps)
				searchTree(pmin, pmax, root->left, cand, eps, trace);
			else
				traceOverL = true;

			if (root->right != NULL && 
				pmin[0] >= root->right->x0 - eps && pmin[1] >= root->right->y0 - eps && pmin[2] >= root->right->z0 - eps && 
				pmin[0] <= root->right->x1 + eps && pmin[1] <= root->right->y1 + eps && pmin[2] <= root->right->z1 + eps &&
				pmax[0] >= root->right->x0 - eps && pmax[1] >= root->right->y0 - eps && pmax[2] >= root->right->z0 - eps && 
				pmax[0] <= root->right->x1 + eps && pmax[1] <= root->right->y1 + eps && pmax[2] <= root->right->z1 + eps)
				searchTree(pmin, pmax, root->right, cand, eps, trace);
			else
				traceOverR = true;
		}

		if ((traceOverL && traceOverR) || trace)
		{
			if (root->left == NULL && root->right == NULL)
			{
				for (int i = 0; i < root->cand.size(); i++)
					cand->push_back(root->cand[i]);
			}
			else
			{
				if (root->left == NULL) searchTree(pmin, pmax, root->right, cand, eps, true);
				if (root->right == NULL) searchTree(pmin, pmax, root->left, cand, eps, true);
			}
		}
	};
	node *construction(vector<float *> range, vector<int> cand)
	{
		if (cand.empty()) return NULL;

		// maximum length
		float r[6];
		for (int i = 0; i < 6; i++) r[i] = range[0][i];

		for (int i = 1; i < range.size(); i++)
		{
			for (int dim = 0; dim < 3; dim++)
			{
				if (r[dim * 2 + 0] > range[i][dim * 2 + 0]) r[dim * 2 + 0] = range[i][dim * 2 + 0];
				if (r[dim * 2 + 1] < range[i][dim * 2 + 1]) r[dim * 2 + 1] = range[i][dim * 2 + 1];
			}
		}

		int dim = 0;
		if (r[3] - r[2] >= r[1] - r[0] && r[3] - r[2] >= r[5] - r[4]) dim = 1;
		if (r[5] - r[4] >= r[1] - r[0] && r[5] - r[4] >= r[3] - r[2]) dim = 2;

		// median
		std::vector<float> data;
		for (int i = 0; i < range.size(); i++)
		{
			data.push_back(range[i][dim * 2 + 0]);
			data.push_back(range[i][dim * 2 + 1]);
		}

		// sort for median
		std::sort(data.begin(), data.end());
		float threshold = data[(int)(data.size() / 2)];

		// classification
		float min, max;
		vector<float *> lrange, rrange;
		vector<int> left, right;
		for (int i = 0; i < range.size(); i++)
		{
			min = range[i][(dim % 3) * 2 + 0];
			max = range[i][(dim % 3) * 2 + 1];
			if (threshold > max)
			{
				left.push_back(cand[i]);
				lrange.push_back(range[i]);
			}
			else if (threshold <= min)
			{
				right.push_back(cand[i]);
				rrange.push_back(range[i]);
			}
			else
			{
				left.push_back(cand[i]);
				right.push_back(cand[i]);
				lrange.push_back(range[i]);
				rrange.push_back(range[i]);
			}
		}

		// new node
		node *elem;
		node *lnode, *rnode;
		if (left.size() == cand.size() || right.size() == cand.size())
		{
			lnode = NULL;
			rnode = NULL;
		}
		else
		{
			lnode = construction(lrange, left);
			rnode = construction(rrange, right);
		}
		if (lnode != NULL && rnode == NULL) elem = lnode;
		else if (lnode == NULL && rnode != NULL) elem = rnode;
		else
		{
			elem = new node();
			elem->left = lnode;
			elem->right = rnode;
			if (lnode == NULL && rnode == NULL) elem->cand = cand;
		}

		return elem;
	};
	void boundingBox(Face f, float *r)
	{
		float minx, maxx, miny, maxy, minz, maxz;
		minx = FLT_MAX, maxx = -FLT_MAX;
		miny = FLT_MAX, maxy = -FLT_MAX;
		minz = FLT_MAX, maxz = -FLT_MAX;
		for (int j = 0; j < 3; j++)
		{
			Vertex v = *f.vertex(j);
			if (minx > v[0]) minx = v[0];
			if (maxx < v[0]) maxx = v[0];
			if (miny > v[1]) miny = v[1];
			if (maxy < v[1]) maxy = v[1];
			if (minz > v[2]) minz = v[2];
			if (maxz < v[2]) maxz = v[2];
		}
		r[0] = minx; r[1] = maxx;
		r[2] = miny; r[3] = maxy;
		r[4] = minz; r[5] = maxz;
	};
	void boundingBox(node *root)
	{
		if (root->left == NULL && root->right == NULL)
		{
			Face f = *m_mesh->face(root->cand[0]);
			float r[6];
			boundingBox(f, r);
			root->x0 = r[0]; root->y0 = r[2]; root->z0 = r[4];
			root->x1 = r[1]; root->y1 = r[3]; root->z1 = r[5];
			for (int i = 1; i < root->cand.size(); i++)
			{
				Face f = *m_mesh->face(root->cand[i]);
				float r[6];
				boundingBox(f, r);
				root->x0 = (root->x0 < r[0]) ? root->x0: r[0];
				root->y0 = (root->y0 < r[2]) ? root->y0: r[2];
				root->z0 = (root->z0 < r[4]) ? root->z0: r[4];
				root->x1 = (root->x1 > r[1]) ? root->x1: r[1];
				root->y1 = (root->y1 > r[3]) ? root->y1: r[3];
				root->z1 = (root->z1 > r[5]) ? root->z1: r[5];
			}
		}
		else if (root->left != NULL && root->right == NULL)
		{
			root->x0 = root->left->x0;
			root->y0 = root->left->y0;
			root->z0 = root->left->z0;
			root->x1 = root->left->x1;
			root->y1 = root->left->y1;
			root->z1 = root->left->z1;
		}
		else if (root->left == NULL && root->right != NULL)
		{
			root->x0 = root->right->x0;
			root->y0 = root->right->y0;
			root->z0 = root->right->z0;
			root->x1 = root->right->x1;
			root->y1 = root->right->y1;
			root->z1 = root->right->z1;
		}
		else
		{
			//assert(root->left != NULL && root->right != NULL);
			root->x0 = (root->left->x0 < root->right->x0) ? root->left->x0: root->right->x0;
			root->y0 = (root->left->y0 < root->right->y0) ? root->left->y0: root->right->y0;
			root->z0 = (root->left->z0 < root->right->z0) ? root->left->z0: root->right->z0;
			root->x1 = (root->left->x1 > root->right->x1) ? root->left->x1: root->right->x1;
			root->y1 = (root->left->y1 > root->right->y1) ? root->left->y1: root->right->y1;
			root->z1 = (root->left->z1 > root->right->z1) ? root->left->z1: root->right->z1;
		}
	};
	void updateTree(node *root)
	{
		if (root->left != NULL) updateTree(root->left);
		if (root->right != NULL) updateTree(root->right);
		
		boundingBox(root);
	};
	node *m_tree;
	Mesh *m_mesh;
};

