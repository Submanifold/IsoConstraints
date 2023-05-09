#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/jet_estimate_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/write_points.h>

#include <utility> // defines std::pair
#include <list>
#include <fstream>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> PointVectorPair;

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

int main(int argc, char*argv[])
{
  const std::string fname = (argc>1) ? argv[1] : CGAL::data_file_path("points_3/sphere_1k.xyz");
  const char* output_fname = (argc > 2) ? argv[2] : "data/fin90_with_PCA_normals_bilateral_smoothed.xyz";
  int nb_neighbors = (argc > 3) ? atoi(argv[3]) : 15; // K-nearest neighbors = 3 rings
  // Reads a point set file in points[].
  std::list<PointVectorPair> points;
  if(!CGAL::IO::read_points(fname,
                            std::back_inserter(points),
                            CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())))
  {
    std::cerr << "Error: cannot read file " << fname<< std::endl;
    return EXIT_FAILURE;
  }
  
  CGAL::jet_estimate_normals<Concurrency_tag>
        (points, nb_neighbors,
         CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
                          .normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));

  if (!CGAL::IO::write_XYZ(output_fname, points,
      CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
      .normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>())
      .stream_precision(8)))
      return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
