#include "ObjIO.h"

#include "iostream"
#include "fstream"
#include "string.h"

using namespace std;

ObjIO::ObjIO(void): MeshIO()
{
}

ObjIO::ObjIO(const char *filename): MeshIO()
{
	read(filename);
}

ObjIO::~ObjIO(void)
{
}

void ObjIO::read(const char *filename)
{
	int nVertex = 0;
	int nFace = 0;
	int nNormal = 0;
	char buf[255];
	
	ifstream fin(filename);
	while (!fin.eof())
	{
		fin.getline(buf, sizeof(buf));
		if (strlen(buf) > 0)
		{
			if (buf[0] == 'v' && buf[1] == ' ') nVertex++;
			if (buf[0] == 'v' && buf[1] == 'n' && buf[2] == ' ') nNormal++;
			if (buf[0] == 'f' && buf[1] == ' ') nFace++;
		}
	}
	if (nNormal != nVertex) nNormal = nVertex;
	else m_hasNormal = true;

	m_nVertex = nVertex;
	m_nFace = nFace;
	m_nNormal = nNormal;

	initArray();

	nVertex = 0;
	nFace = 0;
	nNormal = 0;

	fin.clear();
	fin.seekg(0, ios::beg);
	while (!fin.eof())
	{
		fin.getline(buf, sizeof(buf));
		if (strlen(buf) > 0)
		{
			if (buf[0] == 'v' && buf[1] == ' ')
			{
				sscanf(buf, "%*[^-0-9] %f %f %f", &m_vertex[nVertex * 3], &m_vertex[nVertex * 3 + 1], &m_vertex[nVertex * 3 + 2]);
				nVertex++;
			}
			else if (buf[0] == 'v' && buf[1] == 'n' && buf[2] == ' ')
			{
				sscanf(buf, "%*[^-0-9] %f %f %f", &m_normal[nNormal * 3], &m_normal[nNormal * 3 + 1], &m_normal[nNormal * 3 + 2]);
				nNormal++;
			}
			else if (buf[0] == 'f' && buf[1] == ' ')
			{
				if (strstr(buf, "//") != 0)
				{
					sscanf(buf, "%*[^0-9] %d//%d %d//%d %d//%d", &m_face[nFace * 6], &m_face[nFace * 6 + 3], &m_face[nFace * 6 + 1], &m_face[nFace * 6 + 4], &m_face[nFace * 6 + 2], &m_face[nFace * 6 + 5]);
					m_face[nFace * 6]--;
					m_face[nFace * 6 + 1]--;
					m_face[nFace * 6 + 2]--;
					m_face[nFace * 6 + 3]--;
					m_face[nFace * 6 + 4]--;
					m_face[nFace * 6 + 5]--;
				}
				else
				{
					sscanf(buf, "%*[^0-9] %d %d %d", &m_face[nFace * 3], &m_face[nFace * 3 + 1], &m_face[nFace * 3 + 2]);
					m_face[nFace * 3]--;
					m_face[nFace * 3 + 1]--;
					m_face[nFace * 3 + 2]--;
				}
				nFace++;
			}
		}
	}
	
	fin.close();
}

void ObjIO::save(const char *filename, Mesh *mesh, bool normal)
{
	ofstream fout(filename, ios::out);
	for (int i = 0; i < mesh->nVertex(); i++)
	{
		Vertex v = *mesh->vertex(i);
		fout << "v " << v.fv()[0] << " " << v.fv()[1] << " " << v.fv()[2] << endl;
	}
	if (normal)
	{
		for (int i = 0; i < mesh->nVertex(); i++)
		{
			Normal vn = *mesh->normal(i);
			fout << "vn " << vn.fv()[0] << " " << vn.fv()[1] << " " << vn.fv()[2] << endl;
		}
	}
	for (int i = 0; i < mesh->nFace(); i++)
	{
		Face f = *mesh->face(i);
		if (normal)
			fout << "f " << f.list()[0] + 1 << "//" << f.list()[0] + 1 << " "
						<< f.list()[1] + 1 << "//" << f.list()[1] + 1 << " "
						<< f.list()[2] + 1 << "//" << f.list()[2] + 1 << endl;
		else
			fout << "f " << f.list()[0] + 1 << " " << f.list()[1] + 1 << " " << f.list()[2] + 1<< endl;
	}
	fout.close();
}
