#pragma once
#include "MeshIO.h"

using namespace std;

class ObjIO: public MeshIO
{
public:
	ObjIO(void);
	ObjIO(const char *filename);
	~ObjIO(void);
	void read(const char *filename);
	static void save(const char *filename, Mesh *mesh, bool normal = false);
};