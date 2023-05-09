#include <iostream>
#include <vector>
#include <ctime>
#include <pcl/point_cloud.h>
#include <pcl/octree/octree.h>
#include <boost/thread/thread.hpp>
#include <pcl/visualization/pcl_visualizer.h>
using namespace std;
int main(int argc, char** argv)
{
	if (argc < 4)
	{
		cout << "Require the input, output filename and k value." << endl;
	}
	string input_file_xyz = argv[1];
	string output_index_txt = argv[2];
	string k_string = argv[3];
	int k = stoi(k_string);
	std::ifstream ifs;
	ifs.open(input_file_xyz);
	std::vector<std::vector<double>> points;
	for (int i = 0;; i++)
	{
		vector<double> vec(6, 0); 
		for (auto& d : vec)
		{
			ifs >> d;
		}
		if (ifs.eof())
		{
			break;
		}
		points.push_back(vec);
	}
	std::ofstream ofs;
	ofs.open(output_index_txt, std::ios::ate);

	srand((unsigned int)time(NULL));
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	int num_points = points.size();
	cloud->points.resize(num_points);
	for (size_t i = 0; i < num_points; ++i)
	{
		cloud->points[i].x = points[i][0];
		cloud->points[i].y = points[i][1];
		cloud->points[i].z = points[i][2];
	}
	pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> octree(float(1.0 / 128));
	octree.setInputCloud(cloud);
	octree.addPointsFromInputCloud();
	for (size_t i = 0; i < num_points; ++i)
	{
		pcl::PointXYZ searchPoint;
		searchPoint.x = points[i][0];
		searchPoint.y = points[i][1];
		searchPoint.z = points[i][2];

		vector<int> Idx;
		vector<float> Distance;
		octree.nearestKSearch(searchPoint, k + 1, Idx, Distance);
		
		for (size_t j = 1; j < k; j++) 
		{
			ofs << Idx[j] << " ";
		}
		ofs << Idx[k] << "\n";
	}

	return (0);
}
