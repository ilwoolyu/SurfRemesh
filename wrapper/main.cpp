#include <string>

#include "SRemeshCLP.h"
#include "SurfaceRemeshing.h"

int main(int argc, char* argv[])
{
	PARSE_ARGS;
	
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " --help" << std::endl;
		return -1;
	}
	
	char *subject = (char *) input.c_str();
	char *dfield = NULL;
	if (!coeff.empty()) dfield = (char *) coeff.c_str();
	char *sphere = (char *) temp.c_str();
	char *reference = NULL;
	if (!ref.empty()) reference = (char *) ref.c_str();
	char *dsphere = NULL;
	if (!deform.empty()) dsphere = (char *) deform.c_str();
	char *colormap = NULL;
	if (!color.empty()) colormap = (char *) color.c_str();
	SurfaceRemeshing *SR;

/*	int id = 963368;
	char subject[1024]; sprintf(subject, "/home/hmali/Example/surface/stx_noscale_%d_V12_t1w_label_pp_surf_tMeanSPHARM_procalign.vtk", id);
	char sphere[1024] = "/home/hmali/Example/sphere/stx_noscale_996312_V12_t1w_label_pp_surf_tMeanSPHARM_procalign.sphere.vtk";
	char dfield[1024]; sprintf(dfield, "/NIRAL/work/ilwoolyu/GROUPS-build/result/stx_noscale_%d_V12_t1w_label_pp_surf_tMeanSPHARM_procalign.coeff.txt", id);
	char output[1024]; sprintf(output, "/NIRAL/work/ilwoolyu/SurfRemesh/build/result/stx_noscale_%d_V12_t1w_label_pp_surf_tMeanSPHARM_procalign.vtk", id); */

	cout << "subject: " << subject << endl;
	cout << "sphere: " << sphere << endl;
	if (!coeff.empty())  cout << "deformation: " << dfield << endl;
	cout << "reference: " << reference << endl;
	
	SR = new SurfaceRemeshing(subject, sphere, dfield, keepC, reference, colormap, property);

	cout << "Write output surface model..\n";
	if (!output.empty()) SR->saveDeformedSurface(output.c_str());
	if (dsphere != NULL) SR->saveDeformedSphere(dsphere);
	if (!outputProp.empty()) SR->saveDeformedProperty(outputProp.c_str(), !noheader);

	delete SR;
	
	return 0;
}

