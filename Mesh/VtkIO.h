#pragma once
#include <algorithm>
#include "MeshIO.h"

using namespace std;

class VtkIO: public MeshIO
{
public:
	VtkIO(void);
	VtkIO(const char *filename);
	~VtkIO(void);
	void read(const char *filename);
	static void save(const char *filename, Mesh *mesh, bool normal = false, bool color = false, bool binary = false);
};

