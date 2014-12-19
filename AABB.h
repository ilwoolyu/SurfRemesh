#pragma once
#include <stdio.h>
#include <vector>
#include <algorithm>
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
	int closestFace(float *v, float *coeff)
	{
		float err = 0; // numerical error
		bool success = false;
		int index = -1;

		vector<int> cand;
		float eps = 0;
		for (int trial = 0; cand.empty(); trial++)
		{
			searchTree(v, m_tree, &cand, eps);
			if (trial == 0) eps = 1e-7;
			else eps *= 10;
		}
		sort(cand.begin(), cand.end());
		cand.erase(unique(cand.begin(), cand.end()), cand.end());
		/*for (int i = 0; i < cand.size(); i++)
			cout << cand[i] << " ";
		cout << endl;*/

		float *dist = new float[cand.size()];
		bool *contain = new bool[cand.size()];
		memset(contain, 0, sizeof(bool) * cand.size());

		for (int trial = 0; !success; trial++)
		{
			for (int i = 0; i < cand.size(); i++)
			{
				Face f = *m_mesh->face(cand[i]);
				Vertex a = *f.vertex(0);
				Vertex b = *f.vertex(1);
				Vertex c = *f.vertex(2);

				// bary centric
				Coordinate::cart2bary((float *)a.fv(), (float *)b.fv(), (float *)c.fv(), v, coeff);

				if (coeff[0] >= err && coeff[1] >= err && coeff[2] >= err)
				{
					//printf("%d - (%.10f %.10f %.10f)\n", cand[j] + 1, coeff[0], coeff[1], coeff[2]);
					success = true;
					dist[i] = fabs(Coordinate::dpoint2tri(a.fv(), b.fv(), c.fv(), v));
					contain[i] = true;
				}
			}
			if (trial == 0) err = -1e-16;
			else err *= 10;
		}
		float min_dist = 1;
		for (int i = 0; i < cand.size(); i++)
		{
			if (contain[i] && min_dist > dist[i])
			{
				index = cand[i];
				min_dist = dist[i];
			}
		}
		Face f = *m_mesh->face(index);
		Vertex a = *f.vertex(0);
		Vertex b = *f.vertex(1);
		Vertex c = *f.vertex(2);

		// bary centric
		Coordinate::cart2bary((float *)a.fv(), (float *)b.fv(), (float *)c.fv(), v, coeff);

		delete [] dist;
		delete [] contain;

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
	void searchTree(float *p, node *root, vector<int> *cand, float eps = 0)
	{
		if (root->left != NULL && 
			p[0] >= root->left->x0 - eps && p[1] >= root->left->y0 - eps && p[2] >= root->left->z0 - eps && 
			p[0] <= root->left->x1 + eps && p[1] <= root->left->y1 + eps && p[2] <= root->left->z1 + eps)
			searchTree(p, root->left, cand, eps);
		if (root->right != NULL &&
			p[0] >= root->right->x0 - eps && p[1] >= root->right->y0 - eps && p[2] >= root->right->z0 - eps && 
			p[0] <= root->right->x1 + eps && p[1] <= root->right->y1 + eps && p[2] <= root->right->z1 + eps)
			searchTree(p, root->right, cand, eps);

		if (root->left == NULL && root->right == NULL)
			for (int i = 0; i < root->cand.size(); i++)
				cand->push_back(root->cand[i]);
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

		if (left.size() == cand.size() || right.size() == cand.size()) return NULL;

		// new node
		node *elem;
		node *lnode = construction(lrange, left);
		node *rnode = construction(rrange, right);
		if (lnode != NULL && rnode == NULL) elem = lnode;
		else if (lnode == NULL && rnode != NULL) elem = rnode;
		else
		{
			elem = new node();
			elem->left = lnode;
			elem->right = rnode;
			if (lnode == NULL && rnode == NULL) elem->cand = cand;
		}
		/*node *elem = new node();
		node *lnode = construction(lrange, left, dim + 1);
		node *rnode = construction(rrange, right, dim + 1);
		elem->cand = cand;
		elem->left = lnode;
		elem->right = rnode;*/

		return elem;
	};
	void boundingBox(Face f, float *r)
	{
		float minx, maxx, miny, maxy, minz, maxz;
		minx = 1, maxx = -1;
		miny = 1, maxy = -1;
		minz = 1, maxz = -1;
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

