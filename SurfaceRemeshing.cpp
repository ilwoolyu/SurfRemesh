#include "SurfaceRemeshing.h"
#include "SphericalHarmonics.h"

SurfaceRemeshing::SurfaceRemeshing(void)
{
}

SurfaceRemeshing::~SurfaceRemeshing(void)
{
	delete [] m_x;
	delete [] m_y;
	delete [] m_z;
}

SurfaceRemeshing::SurfaceRemeshing(char *subject, char *sphere, char *dfield, char *sphere_t, vector<string> property, bool interpolation, bool backward)
{
	cout << "loading subject surface model..\n";
	m_subj = new Mesh();
	m_subj->openFile(subject);
	
	int n = m_subj->nVertex();
	m_x = new float[n];
	m_y = new float[n];
	m_z = new float[n];

	for (int i = 0; i < n; i++)
	{
		Vertex *v = (Vertex *)m_subj->vertex(i);
		const float *fv = v->fv();
		m_x[i] = fv[0];
		m_y[i] = fv[1];
		m_z[i] = fv[2];
	}

	// unit sphere information
	cout << "Loading unit sphere information..\n";
	m_sphere_subj = new Mesh();
	m_sphere_subj->openFile(sphere);
	m_sphere_subj->centering();

	if (sphere_t == NULL) sphere_t = sphere;
	m_sphere = new Mesh();
	m_sphere->openFile(sphere_t);
	m_sphere->centering();
	
	m_remesh = new Mesh();
	m_remesh->openFile(sphere_t);

	// forward: subject -> template (color source)
	// backward: template (color source) -> subject
	m_backward = backward;
	m_interpolation = interpolation;
	
	m_property = property;

	FILE *fp;

	if (dfield != NULL)
	{
		cout << "Loading spherical harmonics information..\n";
		// spherical harmonics information
		fp = fopen(dfield, "r");
		fscanf(fp, "%f %f %f", &m_pole[0], &m_pole[1], &m_pole[2]);
		fscanf(fp, "%d", &m_degree);
		int n = (m_degree + 1) * (m_degree + 1);
		m_coeff = new float[n * 2];
		for (int i = 0; i < n; i++)
			fscanf(fp, "%f %f", &m_coeff[i], &m_coeff[n + i]);
		fclose(fp);

		cout << "Initializing deformation..\n";
		// deform vertex
		float *Y = new float[n];
		if (backward)
		{
			for (int i = 0; i < m_sphere->nVertex(); i++)
			{
				float v1[3];
				Vertex *v = (Vertex *)m_sphere->vertex(i);
				SphericalHarmonics::basis(m_degree, (float *)v->fv(), Y);
				reconsCoord(v->fv(), v1, Y, m_coeff, m_degree, m_pole);
				MathVector V(v1); V.unit();
				v->setVertex(V.fv());
			}
		}
		else
		{
			for (int i = 0; i < m_sphere_subj->nVertex(); i++)
			{
				float v1[3];
				Vertex *v = (Vertex *)m_sphere_subj->vertex(i);
				SphericalHarmonics::basis(m_degree, (float *)v->fv(), Y);
				reconsCoord(v->fv(), v1, Y, m_coeff, m_degree, m_pole);
				MathVector V(v1); V.unit();
				v->setVertex(V.fv());
			}
		}
		delete [] Y;
	}
	
	if (property != vector<string>())
	{
		cout << "Load properties..\n";
		int prop = property.size();
		cout << "\t" << prop << " found!\n";
		int nProp = m_sphere_subj->nVertex();
		for (int i = 0; i < prop; i++)
		{
			cout << "\t" << property[i].c_str() << endl;
			fp = fopen(property[i].c_str(), "r");
			char line[1024];
			fgets(line, sizeof(line), fp);
			fgets(line, sizeof(line), fp);
			fgets(line, sizeof(line), fp);
		
			float *refMap = new float[nProp];
			for (int j = 0; j < nProp; j++)
			{
				fscanf(fp, "%f", &refMap[j]);
			}
			fclose(fp);
			m_refMap.push_back(refMap);
		}
	}
	cout << "Tree initialization..\n";
	m_tree = new AABB(m_sphere_subj);

	cout << "Remeshing..\n";
	deformSurface();
	
	if (property != vector<string>())
	{
		cout << "Property transferring..\n";
		deformData();
	}
}

void SurfaceRemeshing::reconsCoord(const float *v0, float *v1, float *Y, float *coeff, float degree, float *pole)
{
	// spharm basis
	int n = (degree + 1) * (degree + 1);

	float mat[9];
	memset(mat, 0, sizeof(float) * 9);
	Coordinate::rotation2equator((float *)v0, mat, pole);

	// rotation to the eqautor
	float rv[3];
	Coordinate::rotPoint(v0, mat, rv);
	
	// polar coodinate
	float phi, theta;
	Coordinate::cart2sph(rv, &phi, &theta);

	// displacement
	float delta[2] = {0, 0};
	for (int i = 0; i < n; i++)
	{
		delta[0] += Y[i] * coeff[i];
		delta[1] += Y[i] * coeff[n + i];
	}
	phi += delta[0]; theta += delta[1];
	Coordinate::sph2cart(phi, theta, rv);

	// inverse rotation
	Coordinate::rotPointInv(rv, mat, v1);
}

float SurfaceRemeshing::dataInterpolation(float *refMap, int index, float *coeff, Mesh *mesh)
{
	float data = 0;
	float err = 0; // numerical error

	Face *f = (Face *)mesh->face(index);
	Vertex *a = (Vertex *)f->vertex(0);
	Vertex *b = (Vertex *)f->vertex(1);
	Vertex *c = (Vertex *)f->vertex(2);
	data = refMap[a->id()] * coeff[0] + refMap[b->id()] * coeff[1] + refMap[c->id()] * coeff[2];

	return data;
}

float SurfaceRemeshing::dataMedian(float *refMap, int index, Mesh *mesh)
{
	float data = 0;

	Face *f = (Face *)mesh->face(index);
	Vertex *a = (Vertex *)f->vertex(0);
	Vertex *b = (Vertex *)f->vertex(1);
	Vertex *c = (Vertex *)f->vertex(2);

	float arrayData[3] = {refMap[a->id()], refMap[b->id()], refMap[c->id()]};
	data = Statistics::median(arrayData, 3);

	return data;
}

void SurfaceRemeshing::deformSurface()
{
	for (int i = 0; i < m_sphere->nVertex(); i++)
	{
		float coeff[3];
		Vertex *v = (Vertex *)m_sphere->vertex(i);
		int id = m_tree->closestFace((float *)v->fv(), coeff);
		float d_v[3];
		if (m_interpolation)
		{
			d_v[0] = dataInterpolation(m_x, id, coeff, m_sphere_subj);
			d_v[1] = dataInterpolation(m_y, id, coeff, m_sphere_subj);
			d_v[2] = dataInterpolation(m_z, id, coeff, m_sphere_subj);
		}
		else
		{
			d_v[0] = dataMedian(m_x, id, m_sphere_subj);
			d_v[1] = dataMedian(m_y, id, m_sphere_subj);
			d_v[2] = dataMedian(m_z, id, m_sphere_subj);
		}
		Vertex *new_v = (Vertex *)m_remesh->vertex(i);
		new_v->setVertex(d_v);
	}
}

void SurfaceRemeshing::deformData()
{
	for (int k = 0; k < m_refMap.size(); k++)
	{
		float *deData = new float[m_sphere->nVertex()];
		for (int i = 0; i < m_sphere->nVertex(); i++)
		{
			float coeff[3];
			Vertex *v = (Vertex *)m_sphere->vertex(i);
			int id = m_tree->closestFace((float *)v->fv(), coeff);
			float d;
			if (m_interpolation)
			{
				d = dataInterpolation(m_refMap[k], id, coeff, m_sphere_subj);
			}
			else
			{
				d = dataMedian(m_refMap[k], id, m_sphere_subj);
			}
			deData[i] = d;
		}
		m_deData.push_back(deData);
	}
}

void SurfaceRemeshing::saveDeformedSurface(char *filename)
{
	m_remesh->saveFile(filename, "vtk");
	if (m_property != vector<string>())
	{
		FILE *fp = fopen(filename, "a");
		fprintf(fp, "FIELD ScalarData %d\n", m_refMap.size());
		for (int i = 0; i < m_refMap.size(); i++)
		{
			string name = m_property[i].substr(0, m_property[i].size() - 4);
			unsigned found = name.find_last_of(".") + 1;
			name = name.substr(found, name.size() - found);
			fprintf(fp, "%s 1 %d float\n", name.c_str(), m_sphere->nVertex());
			for (int j = 0; j < m_sphere->nVertex(); j++)
				fprintf(fp, "%f\n", m_refMap[i][j]);
			fprintf(fp, "\n");
		}
		fclose(fp);
	}
}
