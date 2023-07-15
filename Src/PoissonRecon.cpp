/*
MIT License

Copyright (c) 2023 Dong Xiao, Zuoqiang Shi, Siyu Li, Bailin Deng, Bin Wang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#ifdef _WIN32
#include <Windows.h>
#include <Psapi.h>
#include <list>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <string>
#endif // _WIN32
#include "MyTime.h"
#include "MarchingCubes.h"
#include "Octree.h"
#include "SparseMatrix.h"
#include "CmdLineParser.h"
#include "PPolynomial.h"
#include "Ply.h"
#include "MemoryUsage.h"


#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP
void DumpOutput( const char* format , ... );
void DumpOutput2( char* str , const char* format , ... );
#include "MultiGridOctreeData.h"

#define DEFAULT_FULL_DEPTH 5
#define XSTR(x) STR(x)
#define STR(x) #x
#if DEFAULT_FULL_DEPTH
#pragma message ( "[WARNING] Setting default full depth to " XSTR(DEFAULT_FULL_DEPTH) )
#endif // DEFAULT_FULL_DEPTH

#include <stdarg.h>
char* outputFile=NULL;
int echoStdout=0;
void DumpOutput( const char* format , ... )
{
	if( outputFile )
	{
		FILE* fp = fopen( outputFile , "a" );
		va_list args;
		va_start( args , format );
		vfprintf( fp , format , args );
		fclose( fp );
		va_end( args );
	}
	if( echoStdout )
	{
		va_list args;
		va_start( args , format );
		vprintf( format , args );
		va_end( args );
	}
}
void DumpOutput2( char* str , const char* format , ... )
{
	if( outputFile )
	{
		FILE* fp = fopen( outputFile , "a" );
		va_list args;
		va_start( args , format );
		vfprintf( fp , format , args );
		fclose( fp );
		va_end( args );
	}
	if( echoStdout )
	{
		va_list args;
		va_start( args , format );
		vprintf( format , args );
		va_end( args );
	}
	va_list args;
	va_start( args , format );
	vsprintf( str , format , args );
	va_end( args );
	if( str[strlen(str)-1]=='\n' ) str[strlen(str)-1] = 0;
}


cmdLineString
	In( "in" ) ,
	Out( "out" ) ,
	VoxelGrid( "voxel" ) ,
	XForm( "xForm" );

cmdLineReadable
#ifdef _WIN32
	Performance( "performance" ) ,
#endif // _WIN32
	Complete( "complete") ,
	ShowResidual( "showResidual" ) ,
	NoComments( "noComments" ) ,
	PolygonMesh( "polygonMesh" ) ,
	Confidence( "confidence" ) ,
	NormalWeights( "nWeights" ) ,
	NonManifold( "nonManifold" ) ,
	ASCII( "ascii" ) ,
	Density( "density" ) ,
	Verbose( "verbose" ) ,
	Double( "double" );

cmdLineInt
	Depth( "depth" , 10 ) ,
	CGDepth( "cgDepth" , 0 ) ,
	KernelDepth( "kernelDepth") ,
	AdaptiveExponent( "adaptiveExp" , 1 ) ,
	Iters( "iters" , 8 ) ,
	VoxelDepth( "voxelDepth" , -1 ) ,
	FullDepth( "fullDepth" , DEFAULT_FULL_DEPTH ) ,
	MinDepth( "minDepth" , 0 ) ,
	MaxSolveDepth( "maxSolveDepth" ) ,
	BoundaryType( "boundary" , 1 ) ,
	Threads( "threads" , omp_get_num_procs() );

cmdLineFloat
	SamplesPerNode( "samplesPerNode" , 1.0f ) ,
	Scale( "scale" , 1.1f ) ,
	CSSolverAccuracy( "cgAccuracy" , float(1e-3) ) ,
	PointWeight( "pointWeight" , 10.0f );


cmdLineReadable* params[] =
{
	&In , &Depth , &Out , &XForm ,
	&Scale , &Verbose , &CSSolverAccuracy , &NoComments , &Double ,
	&KernelDepth , &SamplesPerNode , &Confidence , &NormalWeights , &NonManifold , &PolygonMesh , &ASCII , &ShowResidual , &VoxelDepth ,
	&PointWeight , &VoxelGrid , &Threads , &MaxSolveDepth ,
	&AdaptiveExponent , &BoundaryType ,
	&Density ,
	&FullDepth ,
	&MinDepth ,
	&CGDepth , &Iters ,
	&Complete ,
#ifdef _WIN32
	&Performance ,
#endif // _WIN32
};



#ifdef _WIN32
inline double to_seconds( const FILETIME& ft )
{
	const double low_to_sec=100e-9; // 100 nanoseconds
	const double high_to_sec=low_to_sec*4294967296.0;
	return ft.dwLowDateTime*low_to_sec+ft.dwHighDateTime*high_to_sec;
}
#endif // _WIN32

template< class Real, class Vertex >
int PoissonRecon(int depth, float pointweight, std::vector<std::vector<Real>> fine_points, std::vector<std::vector<Real>> fine_normals, std::string file_out, bool oriented_optimized)
{
	std::string oriented_file = "";
	for (int i = 0; i < file_out.length(); i++)
	{
		if (file_out.length() - i >= 5)
		{
			oriented_file = oriented_file + file_out[i];
		}
	}
	if (oriented_optimized)
	{
		oriented_file = oriented_file + "_oriented.xyz";
	}
	else
	{
		oriented_file = oriented_file + "_estimated.xyz";
	}
	std::ofstream ofs_orinted;
	ofs_orinted.open(oriented_file, std::ios::ate);
	for (int i = 0; i < fine_points.size(); i++)
	{
		ofs_orinted << fine_points[i][0] << " " << fine_points[i][1] << " " << fine_points[i][2] << " " << fine_normals[i][0] << " " << fine_normals[i][1] << " " << fine_normals[i][2] << "\n";
	}
	ofs_orinted.close();
	
	int Degree = 2;
	Reset< Real >();

	XForm4x4< Real > xForm, iXForm;
	xForm = XForm4x4< Real >::Identity();
	iXForm = xForm.inverse();

	double t;
	double tt = Time();
	double isoValue = 0;

	Octree< Real > tree;
	tree.threads = Threads.value;
	if (!MaxSolveDepth.set) MaxSolveDepth.value = Depth.value;

	OctNode< TreeNodeData >::SetAllocator(MEMORY_ALLOCATOR_BLOCK_SIZE);

	t = Time();
	int kernelDepth = KernelDepth.set ? KernelDepth.value : Depth.value - 2;
	//int kernelDepth = -1;
	if (kernelDepth > Depth.value)
	{
		fprintf(stderr, "[ERROR] %s can't be greater than %s: %d <= %d\n", KernelDepth.name, Depth.name, KernelDepth.value, Depth.value);
		return EXIT_FAILURE;
	}

	double maxMemoryUsage;
	t = Time(), tree.maxMemoryUsage = 0;
	typename Octree< Real >::PointInfo* pointInfo = new typename Octree< Real >::PointInfo();
	typename Octree< Real >::NormalInfo* normalInfo = new typename Octree< Real >::NormalInfo();
	std::vector< Real >* kernelDensityWeights = new std::vector< Real >();
	std::vector< Real >* centerWeights = new std::vector< Real >();
	int pointCount = tree.template SetTree2< float >(fine_points, fine_normals, MinDepth.value, depth, FullDepth.value, kernelDepth, Real(SamplesPerNode.value), Scale.value, Confidence.set, NormalWeights.set, pointweight, AdaptiveExponent.value, *pointInfo, *normalInfo, *kernelDensityWeights, *centerWeights, BoundaryType.value, xForm, false);
	if (!Density.set) delete kernelDensityWeights, kernelDensityWeights = NULL;

	maxMemoryUsage = tree.maxMemoryUsage;
	t = Time(), tree.maxMemoryUsage = 0;
	Pointer(Real) constraints = tree.SetLaplacianConstraints(*normalInfo);
	delete normalInfo;
	maxMemoryUsage = std::max< double >(maxMemoryUsage, tree.maxMemoryUsage);

	t = Time(), tree.maxMemoryUsage = 0;
	Pointer(Real) solution = tree.SolveSystem(*pointInfo, constraints, ShowResidual.set, Iters.value, MaxSolveDepth.value, CGDepth.value, CSSolverAccuracy.value);
	delete pointInfo;
	FreePointer(constraints);
	maxMemoryUsage = std::max< double >(maxMemoryUsage, tree.maxMemoryUsage);

	CoredFileMeshData< Vertex > mesh;

	if (Verbose.set) tree.maxMemoryUsage = 0;
	t = Time();
	isoValue = tree.GetIsoValue(solution, *centerWeights);
	delete centerWeights;
	tree.GetMCIsoSurface(kernelDensityWeights ? GetPointer(*kernelDensityWeights) : NullPointer< Real >(), solution, isoValue, mesh, true, !NonManifold.set, PolygonMesh.set);
	maxMemoryUsage = std::max< double >(maxMemoryUsage, tree.maxMemoryUsage);
	PlyWritePolygons((char*)file_out.c_str(), &mesh, PLY_ASCII, NULL, 0, iXForm);
	FreePointer(solution);
	return 1;
}


template< class Real, class Vertex >
int ImplicitOrient(int depth, float pointweight, std::string fine_file, std::vector<std::vector<Real>> coarse_points, std::vector<std::vector<Real>> coarse_normals, std::vector<std::vector<Real>> &fine_points, std::vector<std::vector<Real>> &fine_normals)
{
	int Degree = 2;
	Reset< Real >();

	XForm4x4< Real > xForm, iXForm;
	xForm = XForm4x4< Real >::Identity();
	iXForm = xForm.inverse();

	double isoValue = 0;

	Octree< Real > tree;
	tree.threads = Threads.value;
	if (!MaxSolveDepth.set) MaxSolveDepth.value = Depth.value;

	OctNode< TreeNodeData >::SetAllocator(MEMORY_ALLOCATOR_BLOCK_SIZE);

	int kernelDepth = KernelDepth.set ? KernelDepth.value : Depth.value - 2;
	if (kernelDepth > Depth.value)
	{
		fprintf(stderr, "[ERROR] %s can't be greater than %s: %d <= %d\n", KernelDepth.name, Depth.name, KernelDepth.value, Depth.value);
		return EXIT_FAILURE;
	}

	double maxMemoryUsage;
	typename Octree< Real >::PointInfo* pointInfo = new typename Octree< Real >::PointInfo();
	typename Octree< Real >::NormalInfo* normalInfo = new typename Octree< Real >::NormalInfo();
	std::vector< Real >* kernelDensityWeights = new std::vector< Real >();
	std::vector< Real >* centerWeights = new std::vector< Real >();
	int boundarytype = 1;
	int pointCount = tree.template SetTree2< float >(coarse_points, coarse_normals, MinDepth.value, depth, FullDepth.value, kernelDepth, Real(SamplesPerNode.value), Scale.value, Confidence.set, NormalWeights.set, pointweight, AdaptiveExponent.value, *pointInfo, *normalInfo, *kernelDensityWeights, *centerWeights, boundarytype, xForm, false);
	if (!Density.set) delete kernelDensityWeights, kernelDensityWeights = NULL;

	maxMemoryUsage = tree.maxMemoryUsage;
	tree.maxMemoryUsage = 0;
	Pointer(Real) constraints = tree.SetLaplacianConstraints(*normalInfo);

	delete normalInfo;
	maxMemoryUsage = std::max< double >(maxMemoryUsage, tree.maxMemoryUsage);

	tree.maxMemoryUsage = 0;
	Pointer(Real) solution = tree.SolveSystem(*pointInfo, constraints, ShowResidual.set, Iters.value, MaxSolveDepth.value, CGDepth.value, CSSolverAccuracy.value);
	delete pointInfo;
	FreePointer(constraints);
	maxMemoryUsage = std::max< double >(maxMemoryUsage, tree.maxMemoryUsage);

	BSplineData< 2 > _fData;
	_fData.set(tree.tree.maxDepth(), boundarytype);

	std::ifstream ifs_fine;
	ifs_fine.open(fine_file);
	Real _scale = tree._scale;
	Point3D<Real> _center = tree._center;
	std::vector<std::vector<Real>> fine_point_normals;
	while (1)
	{
		std::vector<Real> pn(6, 0);
		ifs_fine >> pn[0] >> pn[1] >> pn[2] >> pn[3] >> pn[4] >> pn[5];
		if (ifs_fine.peek() == EOF)
		{
			break;
		}
		pn[0] = (pn[0] - _center[0]) / _scale;
		pn[1] = (pn[1] - _center[1]) / _scale;
		pn[2] = (pn[2] - _center[2]) / _scale;
		fine_point_normals.push_back(pn);
	}
	int fine_num = fine_point_normals.size();
	std::cout << "The fine point cloud contains " << fine_num << " points." << std::endl;
	if (fine_num == 0)
	{
		std::cout << "File not exist or an empty file." << std::endl;
		exit(1);
	}
	fine_points.resize(fine_num);
	fine_normals.resize(fine_num);
	int threads = omp_get_num_procs();
#pragma omp parallel for num_threads( threads )
	for (int i = 0; i < fine_num; i++)
	{
		fine_points[i].resize(3);
		fine_normals[i].resize(3);
		Real grad_p[3];
		Point3D<Real> p;
		p[0] = fine_point_normals[i][0];
		p[1] = fine_point_normals[i][1];
		p[2] = fine_point_normals[i][2];
		//Calculating the gradient of the implicit field.
		tree.ImplicitGrad(solution, p, grad_p, &(_fData));
		Real dot = grad_p[0] * fine_point_normals[i][3] + grad_p[1] * fine_point_normals[i][4] + grad_p[2] * fine_point_normals[i][5];
		fine_points[i][0] = (p[0] * _scale) + _center[0];
		fine_points[i][1] = (p[1] * _scale) + _center[1];
		fine_points[i][2] = (p[2] * _scale) + _center[2];
		//n*(-1) when setting the tree in SPR, changing back here.
		if (dot >= 0)
		{
			fine_normals[i][0] = -fine_point_normals[i][3];
			fine_normals[i][1] = -fine_point_normals[i][4];
			fine_normals[i][2] = -fine_point_normals[i][5];
		}
		else
		{
			fine_normals[i][0] = fine_point_normals[i][3];
			fine_normals[i][1] = fine_point_normals[i][4];
			fine_normals[i][2] = fine_point_normals[i][5];
		}
	}
	ifs_fine.close();
	FreePointer(solution);
	return 1;
}


//Orienting the coarse point cloud by incorporating isovalue constraints to Poisson equation.
template< class Real, class Vertex >
int Orient_coarse(int depth, std::string file_in, std::vector<std::vector<Real>> &coarse_points, std::vector<std::vector<Real>>& coarse_normals, float alpha, float beta, float pointweight, bool oriented_optimized)
{
	int Degree = 2;
	Reset< Real >();

	XForm4x4< Real > xForm, iXForm;
	xForm = XForm4x4< Real >::Identity();
	iXForm = xForm.inverse();

	double isoValue = 0;
	Octree2< Real > tree;
	tree.threads = Threads.value;
	if (!MaxSolveDepth.set) MaxSolveDepth.value = depth;
	OctNode< TreeNodeData >::SetAllocator(MEMORY_ALLOCATOR_BLOCK_SIZE);
	int kernelDepth = KernelDepth.set ? KernelDepth.value : depth - 2;
	if (kernelDepth > depth)
	{
		fprintf(stderr, "[ERROR] %s can't be greater than %s: %d <= %d\n", KernelDepth.name, Depth.name, KernelDepth.value, Depth.value);
		return EXIT_FAILURE;
	}
	double maxMemoryUsage;
	tree.maxMemoryUsage = 0;
	typename Octree2< Real >::PointInfo* pointInfo = new typename Octree2< Real >::PointInfo();
	typename Octree2< Real >::NormalInfo* normalInfo = new typename Octree2< Real >::NormalInfo();
	typename Octree2< Real >::SmoothCoeff* smoothcoeff = new typename Octree2< Real >::SmoothCoeff();
	std::vector< Real >* kernelDensityWeights = new std::vector< Real >();
	std::vector< Real >* centerWeights = new std::vector< Real >();
	PointStream< float >* pointStream;
	pointStream = new ASCIIPointStream< float >((char*)file_in.c_str());
	int full_depth = 5;
	int pointCount = tree.template SetTree< float >(pointStream, MinDepth.value, depth, full_depth, kernelDepth, Real(SamplesPerNode.value), Scale.value, AdaptiveExponent.value, *pointInfo, *normalInfo, *smoothcoeff, *kernelDensityWeights, *centerWeights, BoundaryType.value, xForm, false);
	if (!Density.set) delete kernelDensityWeights, kernelDensityWeights = NULL;

	maxMemoryUsage = tree.maxMemoryUsage;
	tree.maxMemoryUsage = 0;
	std::vector < std::vector< MatrixEntry <Real> >>  constraints = tree.SetLaplacianConstraints2(*smoothcoeff);
	maxMemoryUsage = std::max< double >(maxMemoryUsage, tree.maxMemoryUsage);
	tree.maxMemoryUsage = 0;
	Pointer(Real) solution = tree.SolveSystem2(file_in, *pointInfo, constraints, coarse_points, coarse_normals, alpha, beta, pointweight, oriented_optimized);
	maxMemoryUsage = std::max< double >(maxMemoryUsage, tree.maxMemoryUsage);
	FreePointer(solution);
}


int main(int argc, char** argv)
{
	//User specific parameters:
	
	//1.coarse_filein: file name of the coarse point cloud containing unoriented normals.
	//2.fine_filein: file name of the fine point cloud containing unoriented normals.
	//3.fileout: file name of the output mesh. The oriented point cloud will be named as fileout[:-4]+"_orinted.xyz".
	//4.use_implicit_orient: whether to use the coarse to fine implicit orientation method. When setting to false, the surface are directly constructed 
	//by the coarse point cloud (after orientation).
	//5.is_noisy_input: whether the input is a noisy point cloud.
	//6.oriented_optimized: "true" for using oriented normals, "false" for using optimized normals (for 3D sketches only)
	

	std::string coarse_filein = ".\\data\\xyzrgb_statuette_coarse.xyz";
	std::string fine_filein = ".\\data\\xyzrgb_statuette_fine.xyz";
	std::string fileout = ".\\results\\xyzrgb_statuette.ply";
	bool use_implicit_orient = true;
	bool is_noisy_input = false;
	bool oriented_optimized = true;

	
	if (argc == 7)
    {
		coarse_filein = argv[1];
		fine_filein = argv[2];
		fileout = argv[3];
        std::string use_implicit_orient_string = argv[4];
		if (use_implicit_orient_string == "true" || use_implicit_orient_string == "True" || use_implicit_orient_string == "TRUE")
		{
			use_implicit_orient = true;
		}
		else if (use_implicit_orient_string == "false" || use_implicit_orient_string == "False" || use_implicit_orient_string == "FALSE")
		{
			use_implicit_orient = false;
		}
		else
		{
			std::cout << "Wrong parameter setting for a bool variable." << std::endl;
			exit(1);
		}
		std::string is_noisy_input_string = argv[5];
		if (is_noisy_input_string == "true" || is_noisy_input_string == "True" || is_noisy_input_string == "TRUE")
		{
			is_noisy_input = true;
		}
		else if (is_noisy_input_string == "false" || is_noisy_input_string == "False" || is_noisy_input_string == "FALSE")
		{
			is_noisy_input = false;
		}
		else
		{
			std::cout << "Wrong parameter setting for a bool variable." << std::endl;
			exit(1);
		}
		std::string oriented_optimized_string = argv[6];
		if (oriented_optimized_string == "true" || oriented_optimized_string == "True" || oriented_optimized_string == "TRUE")
		{
			oriented_optimized = true;
		}
		else if (oriented_optimized_string == "false" || oriented_optimized_string == "False" || oriented_optimized_string == "FALSE")
		{
			oriented_optimized = false;
		}
		else
		{
			std::cout << "Wrong parameter setting for a bool variable." << std::endl;
			exit(1);
		}
    }
	else
	{
		std::cout << "No enough parameters, using the default settings in the main function." << std::endl;
	}
	std::string suffix_coarse = "";
	for (int i = 0; i < coarse_filein.length(); i++)
	{
		if (coarse_filein.length() - i <= 4)
		{
			suffix_coarse = suffix_coarse + coarse_filein[i];
		}
	}
	if (suffix_coarse != ".xyz")
	{
		std::cout << "Only supporting xyz file for coarse point cloud input." << std::endl;
		exit(1);
	}
	if (use_implicit_orient)
	{
		std::string suffix_fine = "";
		for (int i = 0; i < fine_filein.length(); i++)
		{
			if (fine_filein.length() - i <= 4)
			{
				suffix_fine = suffix_fine + fine_filein[i];
			}
		}
		if (suffix_fine != ".xyz")
		{
			std::cout << "Only supporting xyz file for fine point cloud input." << std::endl;
			exit(1);
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	//Numerical parameters, do not need to adjust when running examples.
	float alpha, beta, pointweight;
	if (is_noisy_input)
	{
		alpha = 1e4;
		beta = 1e-2;
		pointweight = 1.0;
	}
	else
	{
		alpha = 1e4;
		beta = 1e-4;
		pointweight = 10.0;
	}
	int coarse_orient_depth = 7;
	int recon_depth = 10;

	std::vector<std::vector<float>> coarse_points;
	std::vector<std::vector<float>> coarse_normals;
	std::vector<std::vector<float>> fine_points;
	std::vector<std::vector<float>> fine_normals;

	//Unoriented normals are required for coarse_filein and fine_filein
	std::ifstream xyz_file;
	xyz_file.open(coarse_filein);
	if (!xyz_file)
	{
		std::cout << "Failed to read coarse xyz file." << std::endl;
		exit(1);
	}

	std::string line;
	while (std::getline(xyz_file, line))
	{
		std::istringstream iss(line);
		float value;
		std::vector<float> values;

		while (iss >> value)
		{
			values.push_back(value);
		}
		if (values.size() != 6)
		{
			std::cout << "Unoriented normals are required for the coarse point cloud, ./jet/normal_estimation.exe can estimate the unoriented normals." << std::endl;
			exit(1);
		}
	}
	xyz_file.close();

	if (use_implicit_orient)
	{
		xyz_file.open(fine_filein);
		if (!xyz_file)
		{
			std::cout << "Failed to read fine xyz file." << std::endl;
			exit(1);
		}

		std::string line;
		while (std::getline(xyz_file, line))
		{
			std::istringstream iss(line);
			float value;
			std::vector<float> values;

			while (iss >> value)
			{
				values.push_back(value);
			}
			if (values.size() != 6)
			{
				std::cout << "Unoriented normals are required for the fine point cloud, ./jet/normal_estimation.exe can estimate the unoriented normals." << std::endl;
				exit(1);
			}
		}
		xyz_file.close();
	}


	std::cout << "Start." << std::endl;
	Orient_coarse<float, PlyValueVertex< float >>(coarse_orient_depth, coarse_filein, coarse_points, coarse_normals, alpha, beta, pointweight, oriented_optimized);
	if (use_implicit_orient)
	{
		std::cout << "Coarse to fine implicit orientation." << std::endl;
		ImplicitOrient<float, PlyValueVertex< float >>(recon_depth, pointweight, fine_filein, coarse_points, coarse_normals, fine_points, fine_normals);
	}
	else
	{
		fine_points = coarse_points;
		fine_normals = coarse_normals;
	}
	PoissonRecon<float, PlyValueVertex< float >>(recon_depth, pointweight, fine_points, fine_normals, fileout, oriented_optimized);
	std::cout << "Finish." << std::endl;
	return 0;
}


