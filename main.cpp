#include "mesh.hpp"
using namespace std;

int main(void)
{
	// reading the high dim mesh with 28k vertices.
	std::string surfMeshFilePath = "/home/sourav/Experimentation/VegaFEM-v3.1/models/sims/male_sitting_stool/mesh.obj";
	

	// std::string surfPointAreaFilePath = "/home/sourav/Experimentation/VegaFEM-v3.1/models/sims/male_sitting_stool/mesh_data/mesh.mr";
	// std::string surfSdfFilePath = "/home/sourav/Experimentation/VegaFEM-v3.1/models/sims/male_sitting_stool/surfmeshviz/sdf_test/mesh.feat";
	std::string  allDeformFilePath = "/home/sourav/Experimentation/VegaFEM-v3.1/models/sims/male_sitting_stool/surfmeshviz/local_area_deform/all_deform.out";
	std::string laplacianFeaturesFilePath = "/home/sourav/Experimentation/VegaFEM-v3.1/models/sims/male_sitting_stool/surfmeshviz/local_area_deform/norm_of_laplacian.feat";
	std::string areaDeformationFilePath = "/home/sourav/Experimentation/VegaFEM-v3.1/models/sims/male_sitting_stool/surfmeshviz/local_area_deform/local_area_deform.feat";
	std::string pointMovementFeatureFilePath = "/home/sourav/Experimentation/VegaFEM-v3.1/models/sims/male_sitting_stool/surfmeshviz/local_area_deform/vertex_movement_dir_norm.feat";

	triangle_mesh_t m(surfMeshFilePath);
	m.loadDeformation(allDeformFilePath);
	// deformation should be specified as space separated input of 3 values.
	// m.laplacianOnDeformations();
	// m.saveAbsLaplacianValuesAsFeature(laplacianFeaturesFilePath);
	// m.saveLaplacianVectorNormAsFeature(laplacianFeaturesFilePath, 4128);
	m.calculateLocalAreaDeformationAboutPoints();
	m.saveLocalAreaDeformationAboutPointsAsFeature(areaDeformationFilePath, 4128);
	m.savePointMovementMagnitudeAsFeature(pointMovementFeatureFilePath);
	// m.saveMassRatio(surfSdfFilePath, surfPointAreaFilePath);
	// m.saveMassRatioAsFeature(surfPointAreaFilePath);

	// std::string deformedMeshFilePath = "/home/sourav/Experimentation/VegaFEM-v3.1/utilities/3dfea/deformedMesh.obj";
	// triangle_mesh_t m_def(deformedMeshFilePath);

	// m - m_def;

	// m.laplacianOnDeformations();
	// m.saveLaplacianVectorNormAsFeature(laplacianFeaturesFilePath, 4128);	

	return 0;
}