#include "VtkIO.h"
#include "fstream"
#include "string.h"
#include "stdlib.h"

using namespace std;

VtkIO::VtkIO(void)
{
}

VtkIO::VtkIO(const char *filename): MeshIO()
{
	read(filename);
}

VtkIO::~VtkIO(void)
{
}

void VtkIO::read(const char *filename)
{
	int nVertex = 0;
	int nFace = 0;
	int nNormal = 0;
	int nColor = 0;
	char buf[16384];
	
	ifstream fin(filename);
	while (!fin.eof())
	{
		fin.getline(buf, sizeof(buf));
		if (strlen(buf) > 0)
		{
			char *tokens;
			char *ptr = strtok(buf, " ");
			if (!strcasecmp(ptr, "points"))
			{
				ptr = strtok(tokens, " ");
				nVertex = atoi(ptr);
				int nElem = 0;
				while (nElem < nVertex * 3)
				{
					fin.getline(buf, sizeof(buf));
					ptr = strtok(buf, " ");
					while (ptr != NULL)
					{
						nElem++;
						ptr = strtok(tokens, " ");
					}
				}
			}
			else if (!strcasecmp(ptr, "polygons"))
			{
				ptr = strtok(tokens, " ");
				nFace = atoi(ptr);
				for (int i = 0; i < nFace; i++)
					fin.getline(buf, sizeof(buf));
			}
			else if (!strcasecmp(ptr, "point_data"))
			{
				ptr = strtok(tokens, " ");
				int n = atoi(ptr);
				fin.getline(buf, sizeof(buf));
				ptr = strtok(buf, " ");
				
				if (!strcasecmp(ptr, "normals"))
				{
					nNormal = n;

					int nElem = 0;
					while (nElem < nNormal * 3)
					{
						fin.getline(buf, sizeof(buf));
						ptr = strtok(buf, " ");
						while (ptr != NULL)
						{
							nElem++;
							ptr = strtok(tokens, " ");
						}
					}
				}
				else if (!strcasecmp(ptr, "color_scalars"))
				{
					m_hasColor = true;
					nColor = n;

					int nElem = 0;
					while (nElem < nColor * 3)
					{
						fin.getline(buf, sizeof(buf));
						ptr = strtok(buf, " ");
						while (ptr != NULL)
						{
							nElem++;
							ptr = strtok(tokens, " ");
						}
					}
				}
			}
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
	nColor = 0;

	fin.clear();
	fin.seekg(0, ios::beg);
	
	while (!fin.eof())
	{
		fin.getline(buf, sizeof(buf));
		if (strlen(buf) > 0)
		{
			char *tokens;
			char *ptr = strtok(buf, " ");
			if (!strcasecmp(ptr, "points"))
			{
				while (nVertex < m_nVertex * 3)
				{
					fin.getline(buf, sizeof(buf));
					ptr = strtok(buf, " ");
					while (ptr != NULL)
					{
						m_vertex[nVertex++] = atof(ptr);
						ptr = strtok(tokens, " ");
					}
				}
			}
			else if (!strcasecmp(ptr, "polygons"))
			{
				for (nFace = 0; nFace < m_nFace; nFace++)
				{
					fin.getline(buf, sizeof(buf));
					if (m_hasNormal)
					{
						sscanf(buf, "%*[0-9] %d %d %d", &m_face[nFace * 6], &m_face[nFace * 6 + 1], &m_face[nFace * 6 + 2]);
						m_face[nFace * 6 + 3] = m_face[nFace * 6];
						m_face[nFace * 6 + 4] = m_face[nFace * 6 + 1];
						m_face[nFace * 6 + 5] = m_face[nFace * 6 + 2];
					}
					else
					{
						sscanf(buf, "%*[0-9] %d %d %d", &m_face[nFace * 3], &m_face[nFace * 3 + 1], &m_face[nFace * 3 + 2]);
					}
				}
			}
			else if (!strcasecmp(ptr, "NORMALS"))
			{
				while (nNormal < m_nNormal * 3)
				{
					fin.getline(buf, sizeof(buf));
					ptr = strtok(buf, " ");
					while (ptr != NULL)
					{
						m_normal[nNormal++] = atof(ptr);
						ptr = strtok(tokens, " ");
					}
				}
			}
			else if (!strcasecmp(ptr, "color_scalars"))
			{
				while (nColor < m_nVertex * 3)
				{
					fin.getline(buf, sizeof(buf));
					ptr = strtok(buf, " ");
					while (ptr != NULL)
					{
						m_color[nColor++] = atof(ptr);
						ptr = strtok(tokens, " ");
					}
				}
			}
		}
	}
	
	fin.close();
}

void VtkIO::save(const char *filename, Mesh *mesh, bool normal, bool color, bool binary)
{
	ofstream fout;
	if (!binary) fout.open(filename, ios::out);
	else fout.open(filename, ios::out | ios::binary);
	fout << "# vtk DataFile Version 3.0" << endl;
	fout << "vtk_output" << endl;
	if (!binary) fout << "ASCII" << endl;
	else fout << "BINARY" << endl;
	fout << "DATASET POLYDATA" << endl;

	fout << "POINTS " << mesh->nVertex() << " float" << endl;
	for (int i = 0; i < mesh->nVertex(); i++)
	{
		const Vertex v = *mesh->vertex(i);
		if (binary)
		{
			float byte[3] = {v.fv()[0], v.fv()[1], v.fv()[2]};
			for (int i = 0; i < 3; i++)
				reverse((char *)&byte[i], (char *)&byte[i] + sizeof(float));
			fout.write((char *)byte, sizeof(float) * 3);
		}
		else
		{
			fout << v.fv()[0] << " " << v.fv()[1] << " " << v.fv()[2];
			fout << endl;
		}
	}
	fout << endl;

	fout << "POLYGONS " << mesh->nFace() << " " << mesh->nFace() * 4 << endl;
	for (int i = 0; i < mesh->nFace(); i++)
	{
		const Face &f = *mesh->face(i);
		if (binary)
		{
			int list = 3;
			reverse((char *)&list, (char *)&list + sizeof(int));
			fout.write((char *)&list, sizeof(int));
			int byte[3] = {f.list()[0], f.list()[1], f.list()[2]};
			for (int j = 0; j < 3; j++)
				reverse((char *)&byte[j], (char *)&byte[j] + sizeof(int));
			fout.write((char *)byte, sizeof(int) * 3);
		}
		else
		{
			fout << "3 " << f.list()[0] << " " << f.list()[1] << " " << f.list()[2] << endl;
		}
	}
	fout << endl;

	if (normal)
	{
		fout << "POINT_DATA " << mesh->nVertex() << endl;
		fout << "NORMALS normals float" << endl;
		for (int i = 0; i < mesh->nVertex(); i++)
		{
			const Vertex &v = *mesh->vertex(i);
			const Normal &vn = *mesh->normal(i);
			if (binary)
			{
				float byte[3] = {vn.fv()[0], vn.fv()[1], vn.fv()[2]};
				for (int j = 0; j < 3; j++)
					reverse((char *)&byte[j], (char *)&byte[j] + sizeof(float));
				fout.write((char *)byte, sizeof(float) * 3);
			}
			else
			{
				Normal vn = *mesh->normal(i);
				fout << vn.fv()[0] << " " << vn.fv()[1] << " " << vn.fv()[2];
				fout << endl;
			}
		}
		fout << endl;
	}
	
	if (color)
	{
		fout << "POINT_DATA " << mesh->nVertex() << endl;
		fout << "COLOR_SCALARS colors 3" << endl;
		for (int i = 0; i < mesh->nVertex(); i++)
		{
			const Vertex &v = *mesh->vertex(i);
			const float *c = v.color();
			if (binary)
			{
				;
			}
			else
			{
				fout << c[0] << " " << c[1] << " " << c[2];
				fout << endl;
			}
		}
		fout << endl;
	}

	fout.close();
}
