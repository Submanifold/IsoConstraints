/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/
#include <iostream>
#include <thread>
#include <windows.h>
#include <cmath>
#include "Octree.h"
#include "MyTime.h"
#include "MemoryUsage.h"
#include "PointStream.h"
#include "MAT.h"
#define ITERATION_POWER 1.0/3
#define MEMORY_ALLOCATOR_BLOCK_SIZE 1<<12
#define SPLAT_ORDER 2

const double MATRIX_ENTRY_EPSILON = 0;
const double EPSILON              = 1e-6;
const double ROUND_EPS            = 1e-5;



//////////////////
// TreeNodeData //
//////////////////
int TreeNodeData::NodeCount = 0;
TreeNodeData::TreeNodeData(void) { nodeIndex = NodeCount++; }
TreeNodeData::~TreeNodeData( void ) { }


////////////
// Octree //
////////////
template< class Real > double Octree< Real >::maxMemoryUsage=0;

template< class Real >
double Octree< Real >::MemoryUsage(void)
{
	double mem = double( MemoryInfo::Usage() ) / (1<<20);
	if( mem>maxMemoryUsage ) maxMemoryUsage=mem;
	return mem;
}

template< class Real >
Octree< Real >::Octree( void )
{
	threads = 1;
	_normalSmooth = 0;
	_constrainValues = false;
}

template< class Real >
bool Octree< Real >::_IsInset( const TreeOctNode* node )
{
	int d , off[3];
	node->depthAndOffset( d , off );
	int res = 1<<d , o = 1<<(d-2);
	return ( off[0]>=o && off[0]<res-o && off[1]>=o && off[1]<res-o && off[2]>=o && off[2]<res-o );
}
template< class Real >
bool Octree< Real >::_IsInsetSupported( const TreeOctNode* node )
{
	int d , off[3];
	node->depthAndOffset( d , off );
	int res = 1<<d , o = (1<<(d-2))-1;
	return ( off[0]>=o && off[0]<res-o && off[1]>=o && off[1]<res-o && off[2]>=o && off[2]<res-o );
}
template< class Real >
int Octree< Real >::SplatOrientedPoint( ConstPointer( Real ) kernelDensityWeights , TreeOctNode* node , const Point3D<Real>& position , const Point3D<Real>& normal , NormalInfo& normalInfo , typename TreeOctNode::NeighborKey3& neighborKey )
{
	double x , dxdy , dxdydz , dx[DIMENSION][SPLAT_ORDER+1];
	double width;
	int off[3];
	typename TreeOctNode::Neighbors3& neighbors = neighborKey.setNeighbors( node );
	Point3D<Real> center;
	Real w;
	node->centerAndWidth( center , w );
	width=w;
	for( int i=0 ; i<3 ; i++ )
	{
#if SPLAT_ORDER==2
		off[i] = 0;
		x = ( center[i] - position[i] - width ) / width;
		dx[i][0] = 1.125+1.500*x+0.500*x*x;
		x = ( center[i] - position[i] ) / width;
		dx[i][1] = 0.750        -      x*x;

		dx[i][2] = 1. - dx[i][1] - dx[i][0];
#elif SPLAT_ORDER==1
		x = ( position[i] - center[i] ) / width;
		if( x<0 )
		{
			off[i] = 0;
			dx[i][0] = -x;
		}
		else
		{
			off[i] = 1;
			dx[i][0] = 1. - x;
		}
		dx[i][1] = 1. - dx[i][0];
#elif SPLAT_ORDER==0
		off[i] = 1;
		dx[i][0] = 1.;
#else
#     error Splat order not supported
#endif // SPLAT_ORDER
	}
	for( int i=off[0] ; i<=off[0]+SPLAT_ORDER ; i++ ) for( int j=off[1] ; j<=off[1]+SPLAT_ORDER ; j++ )
	{
		dxdy = dx[0][i] * dx[1][j];
		for( int k=off[2] ; k<=off[2]+SPLAT_ORDER ; k++ )
			if( neighbors.neighbors[i][j][k] )
			{
				dxdydz = dxdy * dx[2][k];
				TreeOctNode* _node = neighbors.neighbors[i][j][k];
				if( normalInfo.normalIndices.size()<TreeNodeData::NodeCount ) normalInfo.normalIndices.resize( TreeNodeData::NodeCount , -1 );
				int idx = normalInfo.normalIndex( _node );
				if( idx<0 )
				{
					Point3D<Real> n;
					n[0] = n[1] = n[2] = 0;
					idx = normalInfo.normalIndices[ _node->nodeData.nodeIndex ] = (int)normalInfo.normals.size();
					normalInfo.normals.push_back( n );
				}
				normalInfo.normals[idx] += normal * Real( dxdydz );
			}
	}
	return 0;
}
template< class Real >
Real Octree< Real >::SplatOrientedPoint( ConstPointer( Real ) kernelDensityWeights , const Point3D<Real>& position , const Point3D<Real>& normal , NormalInfo& normalInfo , typename TreeOctNode::NeighborKey3& neighborKey , int splatDepth , Real samplesPerNode , int minDepth , int maxDepth )
{
	double dx;
	Point3D<Real> n;
	TreeOctNode* temp;
	int cnt=0;
	double width;
	Point3D< Real > myCenter;
	Real myWidth;
	myCenter[0] = myCenter[1] = myCenter[2] = Real(0.5);
	myWidth = Real(1.0);

	temp = &tree;
	while( temp->depth()<splatDepth )
	{
		if( !temp->children )
		{
			fprintf( stderr , "Octree::SplatOrientedPoint error\n" );
			return -1;
		}
		int cIndex=TreeOctNode::CornerIndex(myCenter,position);
		temp=&temp->children[cIndex];
		myWidth/=2;
		if(cIndex&1) myCenter[0] += myWidth/2;
		else		 myCenter[0] -= myWidth/2;
		if(cIndex&2) myCenter[1] += myWidth/2;
		else		 myCenter[1] -= myWidth/2;
		if(cIndex&4) myCenter[2] += myWidth/2;
		else		 myCenter[2] -= myWidth/2;
	}
	Real weight , depth;
	GetSampleDepthAndWeight( kernelDensityWeights , temp , position , neighborKey , samplesPerNode , depth , weight );

	if( depth<minDepth ) depth=Real(minDepth);
	if( depth>maxDepth ) depth=Real(maxDepth);
	int topDepth=int(ceil(depth));

	dx = 1.0-(topDepth-depth);
	if( topDepth<=minDepth )
	{
		topDepth=minDepth;
		dx=1;
	}
	else if( topDepth>maxDepth )
	{
		topDepth=maxDepth;
		dx=1;
	}
	while( temp->depth()>topDepth ) temp=temp->parent;
	while( temp->depth()<topDepth )
	{
		if(!temp->children) temp->initChildren();
		int cIndex=TreeOctNode::CornerIndex(myCenter,position);
		temp=&temp->children[cIndex];
		myWidth/=2;
		if(cIndex&1) myCenter[0] += myWidth/2;
		else		 myCenter[0] -= myWidth/2;
		if(cIndex&2) myCenter[1] += myWidth/2;
		else		 myCenter[1] -= myWidth/2;
		if(cIndex&4) myCenter[2] += myWidth/2;
		else		 myCenter[2] -= myWidth/2;
	}
	width = 1.0 / ( 1<<temp->depth() );
	n = normal * weight / Real( pow( width , 3 ) ) * Real( dx );
	SplatOrientedPoint( kernelDensityWeights , temp , position , n , normalInfo , neighborKey );
	if( fabs(1.0-dx) > EPSILON )
	{
		dx = Real(1.0-dx);
		temp = temp->parent;
		width = 1.0 / ( 1<<temp->depth() );

		n = normal * weight / Real( pow( width , 3 ) ) * Real( dx );
		SplatOrientedPoint( kernelDensityWeights , temp , position , n , normalInfo , neighborKey );
	}
	return weight;
}

template< class Real >
void Octree< Real >::GetSampleDepthAndWeight( ConstPointer( Real ) kernelDensityWeights , const TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::ConstNeighborKey3& neighborKey , Real samplesPerNode , Real& depth , Real& weight )
{
	const TreeOctNode* temp=node;
	weight = Real(1.0)/GetSampleWeight( kernelDensityWeights , temp , position , neighborKey );
	if( weight>=samplesPerNode ) depth = Real( temp->depth() + log( weight / samplesPerNode ) / log(double(1<<(DIMENSION-1))) );
	else
	{
		Real oldWeight , newWeight;
		oldWeight = newWeight = weight;
		while( newWeight<samplesPerNode && temp->parent )
		{
			temp=temp->parent;
			oldWeight = newWeight;
			newWeight = Real(1.0)/GetSampleWeight( kernelDensityWeights , temp , position, neighborKey );
		}
		depth = Real( temp->depth() + log( newWeight / samplesPerNode ) / log( newWeight / oldWeight ) );
	}
	weight = Real( pow( double(1<<(DIMENSION-1)) , -double(depth) ) );
}
template< class Real >
void Octree< Real >::GetSampleDepthAndWeight( ConstPointer( Real ) kernelDensityWeights , TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::NeighborKey3& neighborKey , Real samplesPerNode , Real& depth , Real& weight )
{
	TreeOctNode* temp=node;
	weight = Real(1.0)/GetSampleWeight( kernelDensityWeights , temp , position , neighborKey );
	if( weight>=samplesPerNode ) depth = Real( temp->depth() + log( weight / samplesPerNode ) / log(double(1<<(DIMENSION-1))) );
	else
	{
		Real oldWeight , newWeight;
		oldWeight = newWeight = weight;
		while( newWeight<samplesPerNode && temp->parent )
		{
			temp=temp->parent;
			oldWeight = newWeight;
			newWeight = Real(1.0)/GetSampleWeight( kernelDensityWeights , temp , position, neighborKey );
		}
		depth = Real( temp->depth() + log( newWeight / samplesPerNode ) / log( newWeight / oldWeight ) );
	}
	weight = Real( pow( double(1<<(DIMENSION-1)) , -double(depth) ) );
}
template< class Real >
Real Octree< Real >::GetSampleWeight( ConstPointer( Real ) kernelDensityWeights , const Point3D<Real>& position , typename TreeOctNode::NeighborKey3& neighborKey , int splatDepth )
{
	Point3D< Real > myCenter;
	Real myWidth;
	myCenter[0] = myCenter[1] = myCenter[2] = Real(0.5);
	myWidth = Real(1.0);

	TreeOctNode* temp = &tree;
	int d = 0;
	while( d<splatDepth )
	{
		if( !temp->children )
		{
			fprintf( stderr , "Octree::SplatOrientedPoint error\n" );
			return -1;
		}
		int cIndex = TreeOctNode::CornerIndex( myCenter , position );
		temp = &temp->children[cIndex];
		myWidth /= 2;
		if( cIndex&1 ) myCenter[0] += myWidth/2;
		else 		   myCenter[0] -= myWidth/2;
		if( cIndex&2 ) myCenter[1] += myWidth/2;
		else 		   myCenter[1] -= myWidth/2;
		if( cIndex&4 ) myCenter[2] += myWidth/2;
		else 		   myCenter[2] -= myWidth/2;
		d++;
	}
	return GetSampleWeight( kernelDensityWeights , temp , position , neighborKey );
}
template< class Real >
Real Octree< Real >::GetSampleWeight( ConstPointer( Real ) kernelDensityWeights , TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::NeighborKey3& neighborKey )
{
	Real weight=0;
	double x , dxdy , dx[DIMENSION][3];
	double width;
	typename TreeOctNode::Neighbors3& neighbors = neighborKey.setNeighbors( node );
	Point3D<Real> center;
	Real w;
	node->centerAndWidth(center,w);
	width=w;

	for( int i=0 ; i<DIMENSION ; i++ )
	{
		x = ( center[i] - position[i] - width ) / width;
		dx[i][0] = 1.125 + 1.500*x + 0.500*x*x;
		x = ( center[i] - position[i] ) / width;
		dx[i][1] = 0.750           -       x*x;

		dx[i][2] = 1.0 - dx[i][1] - dx[i][0];
	}

	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ )
	{
		dxdy = dx[0][i] * dx[1][j];
		for( int k=0 ; k<3 ; k++ ) if( neighbors.neighbors[i][j][k] )
			weight += Real( dxdy * dx[2][k] * kernelDensityWeights[ neighbors.neighbors[i][j][k]->nodeData.nodeIndex ] );
	}
	return Real( 1.0 / weight );
}
template< class Real >
Real Octree< Real >::GetSampleWeight( ConstPointer( Real ) kernelDensityWeights , const TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::ConstNeighborKey3& neighborKey )
{
	Real weight=0;
	double x,dxdy,dx[DIMENSION][3];
	double width;
	typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.getNeighbors( node );
	Point3D<Real> center;
	Real w;
	node->centerAndWidth( center , w );
	width=w;
	//i->dimension j->contribution of two neighbours and itself
	for( int i=0 ; i<DIMENSION ; i++ )
	{
		x = ( center[i] - position[i] - width ) / width;
		dx[i][0] = 1.125 + 1.500*x + 0.500*x*x;
		x = ( center[i] - position[i] ) / width;
		dx[i][1] = 0.750           -       x*x;

		dx[i][2] = 1.0 - dx[i][1] - dx[i][0];
	}

	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ )
	{
		dxdy = dx[0][i] * dx[1][j];
		for( int k=0 ; k<3 ; k++ ) if( neighbors.neighbors[i][j][k] )
			weight += Real( dxdy * dx[2][k] * kernelDensityWeights[ neighbors.neighbors[i][j][k]->nodeData.nodeIndex ] );
	}
	return Real( 1.0 / weight );
}
template< class Real >
int Octree< Real >::UpdateWeightContribution( std::vector< Real >& kernelDensityWeights , TreeOctNode* node , const Point3D<Real>& position , typename TreeOctNode::NeighborKey3& neighborKey , Real weight )
{
	typename TreeOctNode::Neighbors3& neighbors = neighborKey.setNeighbors( node );
	if( kernelDensityWeights.size()<TreeNodeData::NodeCount ) kernelDensityWeights.resize( TreeNodeData::NodeCount , 0 );
	double x , dxdy , dx[DIMENSION][3] , width;
	Point3D< Real > center;
	Real w;
	node->centerAndWidth( center , w );
	width=w;
	const double SAMPLE_SCALE = 1. / ( 0.125 * 0.125 + 0.75 * 0.75 + 0.125 * 0.125 );

	for( int i=0 ; i<DIMENSION ; i++ )
	{
		x = ( center[i] - position[i] - width ) / width;
		dx[i][0] = 1.125 + 1.500*x + 0.500*x*x;
		dx[i][1] = -0.25 - 2.*x - x*x;
		dx[i][2] = 1. - dx[i][1] - dx[i][0];
		// Note that we are splatting along a co-dimension one manifold, so uniform point samples
		// do not generate a unit sample weight.
		dx[i][0] *= SAMPLE_SCALE;
	}
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ )
	{
		dxdy = dx[0][i] * dx[1][j] * weight;
		TreeOctNode** _neighbors = neighbors.neighbors[i][j];
		for( int k=0 ; k<3 ; k++ ) if( _neighbors[k] ) kernelDensityWeights[ _neighbors[k]->nodeData.nodeIndex ] += Real( dxdy * dx[2][k] );
	}
	return 0;
}
template< class Real >
bool Octree< Real >::_InBounds( Point3D< Real > p ) const
{
	if( _boundaryType==0 ){ if( p[0]<Real(0.25) || p[0]>Real(0.75) || p[1]<Real(0.25) || p[1]>Real(0.75) || p[2]<Real(0.25) || p[2]>Real(0.75) ) return false; }
	else                  { if( p[0]<Real(0.00) || p[0]>Real(1.00) || p[1]<Real(0.00) || p[1]>Real(1.00) || p[2]<Real(0.00) || p[2]>Real(1.00) ) return false; }
	return true;
}
template< class Real >
template< class PointReal >
int Octree< Real >::SetTree( PointStream< PointReal >* pointStream , int minDepth , int maxDepth , int fullDepth , 
							int splatDepth , Real samplesPerNode , Real scaleFactor ,
							bool useConfidence , bool useNormalWeights , Real constraintWeight , int adaptiveExponent ,
							PointInfo& pointInfo , NormalInfo& normalInfo , std::vector< Real >& kernelDensityWeights , std::vector< Real >& centerWeights ,
							int boundaryType , XForm4x4< Real > xForm , bool makeComplete )
{
	if( splatDepth<0 ) splatDepth = 0;

	_boundaryType = boundaryType;
	if     ( _boundaryType<0 ) _boundaryType = -1;
	else if( _boundaryType>0 ) _boundaryType =  1;
	_samplesPerNode = samplesPerNode;
	_splatDepth = splatDepth;
	_constrainValues = (constraintWeight>0);

	XForm3x3< Real > xFormN;
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) xFormN(i,j) = xForm(i,j);
	xFormN = xFormN.transpose().inverse();
	minDepth = std::min< int >( minDepth , maxDepth );	// minDepth <= maxDepth
	fullDepth = std::max< int >( minDepth , std::min< int >( fullDepth , maxDepth ) );	// minDepth <= fullDepth <= maxDepth
	// If _boundaryType==0, points are scaled to be in the [0.25,0.75]^3 cube so all depths have to be offset by
	// and the minDepth has to be 2.
	if( _boundaryType==0 )
	{
		minDepth++ , maxDepth++ , fullDepth++;
		if( splatDepth ) splatDepth++;
		minDepth = std::max< int >( minDepth , 2 );
	}
	// Otherwise the points are in the [0,1]^3 cube.
	// However, for Neumann constraints, the function at depth 0 is constant so the system matrix is zero if there
	// is no screening.
#if 0
	else if( _boundaryType==1 && !_constrainValues ) minDepth = std::max< int >( minDepth , 1 );
#endif

	_fData.set( maxDepth , _boundaryType );

	_minDepth = minDepth;
	_fullDepth = fullDepth;
	double pointWeightSum = 0;
	Point3D< Real > min , max , myCenter;
	Real myWidth;
	int i , cnt=0;
	TreeOctNode* temp;

	typename TreeOctNode::NeighborKey3 neighborKey;
	neighborKey.set( maxDepth );

	tree.setFullDepth( _fullDepth );

	// Read through once to get the center and scale
	{
		double t = Time();
		Point3D< Real > p;
		Point3D< PointReal > _p , _n;
		while( pointStream->nextPoint( _p , _n ) )
		{
			p = xForm * Point3D< Real >(_p);
			for( i=0 ; i<DIMENSION ; i++ )
			{
				if( !cnt || p[i]<min[i] ) min[i] = p[i];
				if( !cnt || p[i]>max[i] ) max[i] = p[i];
			}
			cnt++;
		}

		if( _boundaryType==0 ) _scale = std::max< Real >( max[0]-min[0] , std::max< Real >( max[1]-min[1] , max[2]-min[2] ) ) * 2;
		else         _scale = std::max< Real >( max[0]-min[0] , std::max< Real >( max[1]-min[1] , max[2]-min[2] ) );
		_center = ( max+min ) /2;
	}

	_scale *= scaleFactor;
	for( i=0 ; i<DIMENSION ; i++ ) _center[i] -= _scale/2;
	if( splatDepth>0 )
	{
		double t = Time();
		cnt = 0;
		pointStream->reset();
		Point3D< Real > p , n;
		Point3D< PointReal > _p , _n;
		while( pointStream->nextPoint( _p , _n ) )
		{
			p = xForm * Point3D< Real >(_p) , n = xFormN * Point3D< Real >(_n);
			p = ( p - _center ) / _scale;
			if( !_InBounds(p) ) continue;
			myCenter = Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) );
			myWidth = Real(1.0);
			Real weight=Real( 1. );
			if( useConfidence ) weight = Real( Length(n) );
			temp = &tree;
			int d=0;
			while( d<splatDepth )
			{
				//kernelDensityWeight指的是WDˆ (q) 
				UpdateWeightContribution( kernelDensityWeights , temp , p , neighborKey , weight );
				if( !temp->children ) temp->initChildren();
				//cIndex给出p与myCenter的相对位置(x方向大?小?，左前上？右后下？)，判断p在八叉树哪个子节点
				int cIndex=TreeOctNode::CornerIndex( myCenter , p );
				temp = temp->children + cIndex;
				myWidth/=2;
				if( cIndex&1 ) myCenter[0] += myWidth/2;
				else           myCenter[0] -= myWidth/2;
				if( cIndex&2 ) myCenter[1] += myWidth/2;
				else           myCenter[1] -= myWidth/2;
				if( cIndex&4 ) myCenter[2] += myWidth/2;
				else           myCenter[2] -= myWidth/2;
				d++;
			}
			UpdateWeightContribution( kernelDensityWeights , temp , p , neighborKey , weight );
			cnt++;
		}
	}
	kernelDensityWeights.resize( TreeNodeData::NodeCount , 0 );

	std::vector< _PointData >& points = pointInfo.points;

	cnt = 0;
	pointStream->reset();
	Point3D< Real > p , n;
	Point3D< PointReal > _p , _n;
	while( pointStream->nextPoint( _p , _n ) )
	{
		p = xForm * Point3D< Real >(_p) , n = xFormN * Point3D< Real >(_n);
		n *= Real(-1.);
		p = ( p - _center ) / _scale;
		if( !_InBounds(p) ) continue;
		myCenter = Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) );
		myWidth = Real(1.0);
		Real normalLength = Real( Length( n ) );
		if( normalLength!=normalLength || normalLength<=EPSILON ) continue;
		if( !useConfidence ) n /= normalLength;

		Real pointWeight = Real(1.f);
		if( samplesPerNode>0 && splatDepth ) pointWeight = SplatOrientedPoint( GetPointer( kernelDensityWeights ) , p , n , normalInfo , neighborKey , splatDepth , samplesPerNode , _minDepth , maxDepth );
		else
		{
			temp = &tree;
			int d=0;
			if( splatDepth )
			{
				while( d<splatDepth )
				{
					int cIndex=TreeOctNode::CornerIndex(myCenter,p);
					temp = &temp->children[cIndex];
					myWidth /= 2;
					if(cIndex&1) myCenter[0] += myWidth/2;
					else		 myCenter[0] -= myWidth/2;
					if(cIndex&2) myCenter[1] += myWidth/2;
					else		 myCenter[1] -= myWidth/2;
					if(cIndex&4) myCenter[2] += myWidth/2;
					else		 myCenter[2] -= myWidth/2;
					d++;
				}
				pointWeight = GetSampleWeight( GetPointer( kernelDensityWeights ) , temp , p , neighborKey );
			}
			for( i=0 ; i<DIMENSION ; i++ ) n[i] *= pointWeight;
			while( d<maxDepth )
			{
				if( !temp->children ) temp->initChildren();
				int cIndex=TreeOctNode::CornerIndex(myCenter,p);
				temp=&temp->children[cIndex];
				myWidth/=2;
				if(cIndex&1) myCenter[0] += myWidth/2;
				else		 myCenter[0] -= myWidth/2;
				if(cIndex&2) myCenter[1] += myWidth/2;
				else		 myCenter[1] -= myWidth/2;
				if(cIndex&4) myCenter[2] += myWidth/2;
				else		 myCenter[2] -= myWidth/2;
				d++;
			}
			SplatOrientedPoint( GetPointer( kernelDensityWeights ) , temp , p , n , normalInfo , neighborKey );
		}
		pointWeightSum += pointWeight;
		//_constrainValues是bool变量，表明是否使用SPR中的point constrains
		if( _constrainValues )
		{
			Real pointScreeningWeight = useNormalWeights ? Real( normalLength ) : Real(1.f);
			int d = 0;
			TreeOctNode* temp = &tree;
			myCenter = Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) );
			myWidth = Real(1.0);
			while( 1 )
			{
				if( pointInfo.pointIndices.size()<TreeNodeData::NodeCount ) pointInfo.pointIndices.resize( TreeNodeData::NodeCount , -1 );
				int idx = pointInfo.pointIndex( temp );

				if( idx==-1 )
				{
					idx = (int)points.size();
					points.push_back( _PointData( p*pointScreeningWeight , pointScreeningWeight ) );
					pointInfo.pointIndices[ temp->nodeData.nodeIndex ] = idx;
				}
				else
				{
					points[idx].weight += pointScreeningWeight;
					points[idx].position += p*pointScreeningWeight;
				}

				int cIndex = TreeOctNode::CornerIndex( myCenter , p );
				if( !temp->children ) break;
				temp = &temp->children[cIndex];
				myWidth /= 2;
				if( cIndex&1 ) myCenter[0] += myWidth/2;
				else		   myCenter[0] -= myWidth/2;
				if( cIndex&2 ) myCenter[1] += myWidth/2;
				else		   myCenter[1] -= myWidth/2;
				if( cIndex&4 ) myCenter[2] += myWidth/2;
				else		   myCenter[2] -= myWidth/2;
				d++;
			}
		}
		cnt++;
	}

	if( _boundaryType==0 ) pointWeightSum *= Real(4.);
	constraintWeight *= Real( pointWeightSum );
	constraintWeight /= cnt;

	MemoryUsage( );
	if( _constrainValues )
		// Set the average position and scale the weights
		for( TreeOctNode* node=tree.nextNode() ; node ; node=tree.nextNode(node) )
			if( pointInfo.pointIndex( node )!=-1 )
			{
				int idx = pointInfo.pointIndex( node );
				points[idx].position /= points[idx].weight;
				int e = ( _boundaryType==0 ? node->depth()-1 : node->depth() ) * adaptiveExponent - ( _boundaryType==0 ? maxDepth-1 : maxDepth ) * (adaptiveExponent-1);
				if( e<0 ) points[idx].weight /= Real( 1<<(-e) );
				else      points[idx].weight *= Real( 1<<  e  );
				points[idx].weight *= Real( constraintWeight );
			}
#if FORCE_NEUMANN_FIELD
	if( _boundaryType==1 )
		for( TreeOctNode* node=tree.nextNode() ; node ; node=tree.nextNode( node ) )
		{
			int d , off[3] , res;
			node->depthAndOffset( d , off );
			res = 1<<d;
			int idx = normalInfo.normalIndex( node );
			if( idx<0 ) continue;
			Point3D< Real >& normal = normalInfo.normals[ idx ];
			for( int d=0 ; d<3 ; d++ ) if( off[d]==0 || off[d]==res-1 ) normal[d] = 0;
		}
#endif // FORCE_NEUMANN_FIELD
	centerWeights.resize( tree.nodes() , 0 );
	kernelDensityWeights.resize( tree.nodes() , 0 );
	// Set the point weights for evaluating the iso-value
	//centerWeights = Length( normalInfo.normals[ idx ] )记为W_D(s.p)
	for( TreeOctNode* node=tree.nextNode() ; node ; node=tree.nextNode(node) )
	{
		int idx = normalInfo.normalIndex( node );
		if( idx<0 ) centerWeights[ node->nodeData.nodeIndex ] = 0;
		else        centerWeights[ node->nodeData.nodeIndex ] = Real( Length( normalInfo.normals[ idx ] ) );
	}
	MemoryUsage();
	{
		std::vector< int > indexMap;
		if( makeComplete ) MakeComplete( &indexMap );
		else ClipTree( normalInfo ) , Finalize( &indexMap );

		{
			std::vector< int > temp = pointInfo.pointIndices;
			pointInfo.pointIndices.resize( indexMap.size() );
			for( int i=0 ; i<indexMap.size() ; i++ )
				if( indexMap[i]<temp.size() ) pointInfo.pointIndices[i] = temp[ indexMap[i] ];
				else                          pointInfo.pointIndices[i] = -1;
		}
		{
			std::vector< int > temp = normalInfo.normalIndices;
			normalInfo.normalIndices.resize( indexMap.size() );
			for( int i=0 ; i<indexMap.size() ; i++ )
				if( indexMap[i]<temp.size() ) normalInfo.normalIndices[i] = temp[ indexMap[i] ];
				else                          normalInfo.normalIndices[i] = -1;
		}
		{
			std::vector< Real > temp = centerWeights;
			centerWeights.resize( indexMap.size() );
			for( int i=0 ; i<indexMap.size() ; i++ )
				if( indexMap[i]<temp.size() ) centerWeights[i] = temp[ indexMap[i] ];
				else                          centerWeights[i] = (Real)0;
		}
		{
			std::vector< Real > temp = kernelDensityWeights;
			kernelDensityWeights.resize( indexMap.size() );
			for( int i=0 ; i<indexMap.size() ; i++ )
				if( indexMap[i]<temp.size() ) kernelDensityWeights[i] = temp[ indexMap[i] ];
				else                          kernelDensityWeights[i] = (Real)0;
		}
	}
	return cnt;
}
template< class Real >
template< class PointReal >
int Octree< Real >::SetTree2(std::vector<std::vector<Real>> in_points, std::vector<std::vector<Real>> in_normals, int minDepth, int maxDepth, int fullDepth,
	int splatDepth, Real samplesPerNode, Real scaleFactor,
	bool useConfidence, bool useNormalWeights, Real constraintWeight, int adaptiveExponent,
	PointInfo& pointInfo, NormalInfo& normalInfo, std::vector< Real >& kernelDensityWeights, std::vector< Real >& centerWeights,
	int boundaryType, XForm4x4< Real > xForm, bool makeComplete)
{
	if (splatDepth < 0) splatDepth = 0;

	_boundaryType = boundaryType;
	if (_boundaryType < 0) _boundaryType = -1;
	else if (_boundaryType > 0) _boundaryType = 1;
	_samplesPerNode = samplesPerNode;
	_splatDepth = splatDepth;
	_constrainValues = (constraintWeight > 0);
	XForm3x3< Real > xFormN;
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) xFormN(i, j) = xForm(i, j);
	xFormN = xFormN.transpose().inverse();
	minDepth = std::min< int >(minDepth, maxDepth);	// minDepth <= maxDepth
	fullDepth = std::max< int >(minDepth, std::min< int >(fullDepth, maxDepth));	// minDepth <= fullDepth <= maxDepth
	// If _boundaryType==0, points are scaled to be in the [0.25,0.75]^3 cube so all depths have to be offset by
	// and the minDepth has to be 2.
	if (_boundaryType == 0)
	{
		minDepth++, maxDepth++, fullDepth++;
		if (splatDepth) splatDepth++;
		minDepth = std::max< int >(minDepth, 2);
	}
	// Otherwise the points are in the [0,1]^3 cube.
	// However, for Neumann constraints, the function at depth 0 is constant so the system matrix is zero if there
	// is no screening.
#if 0
	else if (_boundaryType == 1 && !_constrainValues) minDepth = std::max< int >(minDepth, 1);
#endif

	_fData.set(maxDepth, _boundaryType);

	_minDepth = minDepth;
	_fullDepth = fullDepth;
	double pointWeightSum = 0;
	Point3D< Real > min, max, myCenter;
	Real myWidth;
	int i, cnt = 0;
	TreeOctNode* temp;

	typename TreeOctNode::NeighborKey3 neighborKey;
	neighborKey.set(maxDepth);

	tree.setFullDepth(_fullDepth);

	// Read through once to get the center and scale
	{
		double t = Time();
		Point3D< Real > p;
		for (int ii = 0; ii < in_points.size(); ii++)
		{
			Point3D< Real > _p;
			_p[0] = in_points[ii][0];
			_p[1] = in_points[ii][1];
			_p[2] = in_points[ii][2];
			p = xForm * Point3D< Real >(_p);
			for (i = 0; i < DIMENSION; i++)
			{
				if (!cnt || p[i] < min[i]) min[i] = p[i];
				if (!cnt || p[i] > max[i]) max[i] = p[i];
			}
			cnt++;
		}

		if (_boundaryType == 0) _scale = std::max< Real >(max[0] - min[0], std::max< Real >(max[1] - min[1], max[2] - min[2])) * 2;
		else         _scale = std::max< Real >(max[0] - min[0], std::max< Real >(max[1] - min[1], max[2] - min[2]));
		_center = (max + min) / 2;
	}
	_scale *= scaleFactor;
	for (i = 0; i < DIMENSION; i++) _center[i] -= _scale / 2;
	if (splatDepth > 0)
	{
		double t = Time();
		cnt = 0;

		Point3D< Real > p, n;

		for (int ii = 0; ii < in_points.size(); ii++)
		{
			Point3D< Real > _p, _n;
			_p[0] = in_points[ii][0];
			_p[1] = in_points[ii][1];
			_p[2] = in_points[ii][2];
			_n[0] = in_normals[ii][0];
			_n[1] = in_normals[ii][1];
			_n[2] = in_normals[ii][2];
			p = xForm * Point3D< Real >(_p), n = xFormN * Point3D< Real >(_n);
			p = (p - _center) / _scale;
			if (!_InBounds(p)) continue;
			myCenter = Point3D< Real >(Real(0.5), Real(0.5), Real(0.5));
			myWidth = Real(1.0);
			Real weight = Real(1.);
			if (useConfidence) weight = Real(Length(n));
			temp = &tree;
			int d = 0;
			while (d < splatDepth)
			{
				UpdateWeightContribution(kernelDensityWeights, temp, p, neighborKey, weight);
				if (!temp->children) temp->initChildren();
				int cIndex = TreeOctNode::CornerIndex(myCenter, p);
				temp = temp->children + cIndex;
				myWidth /= 2;
				if (cIndex & 1) myCenter[0] += myWidth / 2;
				else           myCenter[0] -= myWidth / 2;
				if (cIndex & 2) myCenter[1] += myWidth / 2;
				else           myCenter[1] -= myWidth / 2;
				if (cIndex & 4) myCenter[2] += myWidth / 2;
				else           myCenter[2] -= myWidth / 2;
				d++;
			}
			UpdateWeightContribution(kernelDensityWeights, temp, p, neighborKey, weight);
			cnt++;
		}
	}
	kernelDensityWeights.resize(TreeNodeData::NodeCount, 0);

	std::vector< _PointData >& points = pointInfo.points;

	cnt = 0;
	Point3D< Real > p, n;
	for (int ii = 0; ii < in_points.size(); ii++)
	{
		Point3D< Real > _p, _n;
		_p[0] = in_points[ii][0];
		_p[1] = in_points[ii][1];
		_p[2] = in_points[ii][2];
		_n[0] = in_normals[ii][0];
		_n[1] = in_normals[ii][1];
		_n[2] = in_normals[ii][2];
		p = xForm * Point3D< Real >(_p), n = xFormN * Point3D< Real >(_n);
		n *= Real(-1.);
		p = (p - _center) / _scale;

		if (!_InBounds(p)) continue;
		myCenter = Point3D< Real >(Real(0.5), Real(0.5), Real(0.5));
		myWidth = Real(1.0);
		Real normalLength = Real(Length(n));
		if (normalLength != normalLength || normalLength <= EPSILON) continue;
		if (!useConfidence) n /= normalLength;

		Real pointWeight = Real(1.f);
		if (samplesPerNode > 0 && splatDepth) pointWeight = SplatOrientedPoint(GetPointer(kernelDensityWeights), p, n, normalInfo, neighborKey, splatDepth, samplesPerNode, _minDepth, maxDepth);
		else
		{
			temp = &tree;
			int d = 0;
			if (splatDepth)
			{
				while (d < splatDepth)
				{
					int cIndex = TreeOctNode::CornerIndex(myCenter, p);
					temp = &temp->children[cIndex];
					myWidth /= 2;
					if (cIndex & 1) myCenter[0] += myWidth / 2;
					else		 myCenter[0] -= myWidth / 2;
					if (cIndex & 2) myCenter[1] += myWidth / 2;
					else		 myCenter[1] -= myWidth / 2;
					if (cIndex & 4) myCenter[2] += myWidth / 2;
					else		 myCenter[2] -= myWidth / 2;
					d++;
				}
				pointWeight = GetSampleWeight(GetPointer(kernelDensityWeights), temp, p, neighborKey);
			}
			for (i = 0; i < DIMENSION; i++) n[i] *= pointWeight;
			while (d < maxDepth)
			{
				if (!temp->children) temp->initChildren();
				int cIndex = TreeOctNode::CornerIndex(myCenter, p);
				temp = &temp->children[cIndex];
				myWidth /= 2;
				if (cIndex & 1) myCenter[0] += myWidth / 2;
				else		 myCenter[0] -= myWidth / 2;
				if (cIndex & 2) myCenter[1] += myWidth / 2;
				else		 myCenter[1] -= myWidth / 2;
				if (cIndex & 4) myCenter[2] += myWidth / 2;
				else		 myCenter[2] -= myWidth / 2;
				d++;
			}
			SplatOrientedPoint(GetPointer(kernelDensityWeights), temp, p, n, normalInfo, neighborKey);
		}
		pointWeightSum += pointWeight;
		if (_constrainValues)
		{
			Real pointScreeningWeight = useNormalWeights ? Real(normalLength) : Real(1.f);
			int d = 0;
			TreeOctNode* temp = &tree;
			myCenter = Point3D< Real >(Real(0.5), Real(0.5), Real(0.5));
			myWidth = Real(1.0);
			while (1)
			{
				if (pointInfo.pointIndices.size() < TreeNodeData::NodeCount) pointInfo.pointIndices.resize(TreeNodeData::NodeCount, -1);
				int idx = pointInfo.pointIndex(temp);

				if (idx == -1)
				{
					idx = (int)points.size();
					points.push_back(_PointData(p * pointScreeningWeight, pointScreeningWeight));
					pointInfo.pointIndices[temp->nodeData.nodeIndex] = idx;
				}
				else
				{
					points[idx].weight += pointScreeningWeight;
					points[idx].position += p * pointScreeningWeight;
				}

				int cIndex = TreeOctNode::CornerIndex(myCenter, p);
				if (!temp->children) break;
				temp = &temp->children[cIndex];
				myWidth /= 2;
				if (cIndex & 1) myCenter[0] += myWidth / 2;
				else		   myCenter[0] -= myWidth / 2;
				if (cIndex & 2) myCenter[1] += myWidth / 2;
				else		   myCenter[1] -= myWidth / 2;
				if (cIndex & 4) myCenter[2] += myWidth / 2;
				else		   myCenter[2] -= myWidth / 2;
				d++;
			}
		}
		cnt++;
	}
	if (_boundaryType == 0) pointWeightSum *= Real(4.);
	constraintWeight *= Real(pointWeightSum);
	constraintWeight /= cnt;

	MemoryUsage();
	if (_constrainValues)
		// Set the average position and scale the weights
		for (TreeOctNode* node = tree.nextNode(); node; node = tree.nextNode(node))
			if (pointInfo.pointIndex(node) != -1)
			{
				int idx = pointInfo.pointIndex(node);
				points[idx].position /= points[idx].weight;
				int e = (_boundaryType == 0 ? node->depth() - 1 : node->depth()) * adaptiveExponent - (_boundaryType == 0 ? maxDepth - 1 : maxDepth) * (adaptiveExponent - 1);
				if (e < 0) points[idx].weight /= Real(1 << (-e));
				else      points[idx].weight *= Real(1 << e);
				points[idx].weight *= Real(constraintWeight);
			}
#if FORCE_NEUMANN_FIELD
	if (_boundaryType == 1)
		for (TreeOctNode* node = tree.nextNode(); node; node = tree.nextNode(node))
		{
			int d, off[3], res;
			node->depthAndOffset(d, off);
			res = 1 << d;
			int idx = normalInfo.normalIndex(node);
			if (idx < 0) continue;
			Point3D< Real >& normal = normalInfo.normals[idx];
			for (int d = 0; d < 3; d++)
			{
				if (off[d] == 0 || off[d] == res - 1)
				{
					normal[d] = 0;
				}
			}

		}
#endif // FORCE_NEUMANN_FIELD
	centerWeights.resize(tree.nodes(), 0);
	kernelDensityWeights.resize(tree.nodes(), 0);
	//Set the point weights for evaluating the iso-value
	//node处的normalInfo不为空，则centerWeights是length，否则为0
	for (TreeOctNode* node = tree.nextNode(); node; node = tree.nextNode(node))
	{
		int idx = normalInfo.normalIndex(node);
		if (idx < 0) centerWeights[node->nodeData.nodeIndex] = 0;
		else        centerWeights[node->nodeData.nodeIndex] = Real(Length(normalInfo.normals[idx]));
	}
	MemoryUsage();
	{
		std::vector< int > indexMap;
		if (makeComplete) MakeComplete(&indexMap);
		else ClipTree(normalInfo), Finalize(&indexMap);

		{
			std::vector< int > temp = pointInfo.pointIndices;
			pointInfo.pointIndices.resize(indexMap.size());
			for (int i = 0; i < indexMap.size(); i++)
				if (indexMap[i] < temp.size()) pointInfo.pointIndices[i] = temp[indexMap[i]];
				else                          pointInfo.pointIndices[i] = -1;
		}
		{
			std::vector< int > temp = normalInfo.normalIndices;
			normalInfo.normalIndices.resize(indexMap.size());
			for (int i = 0; i < indexMap.size(); i++)
				if (indexMap[i] < temp.size()) normalInfo.normalIndices[i] = temp[indexMap[i]];
				else                          normalInfo.normalIndices[i] = -1;
		}
		{
			std::vector< Real > temp = centerWeights;
			centerWeights.resize(indexMap.size());
			for (int i = 0; i < indexMap.size(); i++)
				if (indexMap[i] < temp.size()) centerWeights[i] = temp[indexMap[i]];
				else                          centerWeights[i] = (Real)0;
		}
		{
			std::vector< Real > temp = kernelDensityWeights;
			kernelDensityWeights.resize(indexMap.size());
			for (int i = 0; i < indexMap.size(); i++)
				if (indexMap[i] < temp.size()) kernelDensityWeights[i] = temp[indexMap[i]];
				else                          kernelDensityWeights[i] = (Real)0;
		}
	}
	return cnt;
}

template< class Real >
void Octree< Real >::MakeComplete( std::vector< int >* map )
{
	tree.setFullDepth( tree.maxDepth() );
	refineBoundary( map );
	MemoryUsage();
}
template< class Real >
void Octree< Real >::ClipTree( const NormalInfo& normalInfo )
{
	int maxDepth = tree.maxDepth();
	for( TreeOctNode* temp=tree.nextNode() ; temp ; temp=tree.nextNode(temp) )
		if( temp->children && temp->depth()>=_fullDepth )
		{
			int hasNormals=0;
			for( int i=0 ; i<Cube::CORNERS && !hasNormals ; i++ ) hasNormals = HasNormals( &temp->children[i] , normalInfo );
			if( !hasNormals ) temp->children=NULL;
		}
	MemoryUsage();
}

template< class Real >
void Octree< Real >::Finalize( std::vector< int >* map )
{
	int maxDepth = tree.maxDepth( );
	typename TreeOctNode::NeighborKey3 neighborKey;
	neighborKey.set( maxDepth );
	for( int d=maxDepth ; d>1 ; d-- )
		for( TreeOctNode* node=tree.nextNode() ; node ; node=tree.nextNode( node ) ) if( node->depth()==d )
		{
			typename TreeOctNode::Neighbors3& neighbors = neighborKey.setNeighbors( node->parent->parent );
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
				if( neighbors.neighbors[i][j][k] && !neighbors.neighbors[i][j][k]->children )
					neighbors.neighbors[i][j][k]->initChildren();
		}
	refineBoundary( map );
}
template< class Real >
double Octree< Real >::GetLaplacian( const typename BSplineData< 2 >::Integrator& integrator , int d , const int off1[] , const int off2[] , bool childParent ) const
{
	double vv[] =
	{
		integrator.dot( d , off1[0] , off2[0] , false , false , childParent ) ,
		integrator.dot( d , off1[1] , off2[1] , false , false , childParent ) ,
		integrator.dot( d , off1[2] , off2[2] , false , false , childParent )
	};
	double dd[] =
	{
		integrator.dot( d , off1[0] , off2[0] , true , true , childParent ) ,
		integrator.dot( d , off1[1] , off2[1] , true , true , childParent ) ,
		integrator.dot( d , off1[2] , off2[2] , true , true , childParent )
	};
	return dd[0]*vv[1]*vv[2] + vv[0]*dd[1]*vv[2] + vv[0]*vv[1]*dd[2];
}
template< class Real >
double Octree< Real >::GetDivergence1( const typename BSplineData< 2 >::Integrator& integrator , int d , const int off1[] , const int off2[] , bool childParent , const Point3D< Real >& normal1 ) const
{
	return Point3D< double >::Dot( GetDivergence1( integrator , d , off1 , off2 , childParent ) , normal1 );
}
template< class Real > 
double Octree< Real >::GetDivergence2( const typename BSplineData< 2 >::Integrator& integrator , int d , const int off1[] , const int off2[] , bool childParent , const Point3D< Real >& normal2 ) const
{
	return Point3D< double >::Dot( GetDivergence2( integrator , d , off1 , off2 , childParent ) , normal2 );
}
template< class Real >
Point3D< double > Octree< Real >::GetDivergence1( const typename BSplineData< 2 >::Integrator& integrator , int d , const int off1[] , const int off2[] , bool childParent ) const
{
	//B(p) = B(x)B(y)B(z), grad(B(p)) = grad(B(x))B(y)B(z) + ... + ...
	double vv[] =
	{
		integrator.dot( d , off1[0] , off2[0] , false , false , childParent ) ,
		integrator.dot( d , off1[1] , off2[1] , false , false , childParent ) ,
		integrator.dot( d , off1[2] , off2[2] , false , false , childParent )
	};
#if GRADIENT_DOMAIN_SOLUTION
	// Take the dot-product of the vector-field with the gradient of the basis function
	double vd[] = 
	{
		integrator.dot( d , off1[0] , off2[0] , false , true , childParent ) ,
		integrator.dot( d , off1[1] , off2[1] , false , true , childParent ) ,
		integrator.dot( d , off1[2] , off2[2] , false , true , childParent )
	};
	return  Point3D< double >( vd[0]*vv[1]*vv[2] , vv[0]*vd[1]*vv[2] , vv[0]*vv[1]*vd[2] );
#else // !GRADIENT_DOMAIN_SOLUTION
	// Take the dot-product of the divergence of the vector-field with the basis function
	double dv[] = 
	{
		integrator.dot( d , off1[0] , off2[0] , true , false , childParent ) ,
		integrator.dot( d , off1[1] , off2[1] , true , false , childParent ) ,
		integrator.dot( d , off1[2] , off2[2] , true , false , childParent )
	};
	return  -Point3D< double >( dv[0]*vv[1]*vv[2] , vv[0]*dv[1]*vv[2] , vv[0]*vv[1]*dv[2] );
#endif // GRADIENT_DOMAIN_SOLUTION
}
template< class Real > 
Point3D< double > Octree< Real >::GetDivergence2( const typename BSplineData< 2 >::Integrator& integrator , int d , const int off1[] , const int off2[] , bool childParent ) const
{
	double vv[] =
	{
		integrator.dot( d , off1[0] , off2[0] , false , false , childParent ) ,
		integrator.dot( d , off1[1] , off2[1] , false , false , childParent ) ,
		integrator.dot( d , off1[2] , off2[2] , false , false , childParent )
	};
#if GRADIENT_DOMAIN_SOLUTION
	// Take the dot-product of the vector-field with the gradient of the basis function
	double dv[] = 
	{
		integrator.dot( d , off1[0] , off2[0] , true , false , childParent ) ,
		integrator.dot( d , off1[1] , off2[1] , true , false , childParent ) ,
		integrator.dot( d , off1[2] , off2[2] , true , false , childParent )
	};
	return  Point3D< double >( dv[0]*vv[1]*vv[2] , vv[0]*dv[1]*vv[2] , vv[0]*vv[1]*dv[2] );
#else // !GRADIENT_DOMAIN_SOLUTION
	// Take the dot-product of the divergence of the vector-field with the basis function
	double vd[] = 
	{
		integrator.dot( d , off1[0] , off2[0] , false , true , childParent ) ,
		integrator.dot( d , off1[1] , off2[1] , false , true , childParent ) ,
		integrator.dot( d , off1[2] , off2[2] , false , true , childParent )
	};
	return -Point3D< double >( vd[0]*vv[1]*vv[2] , vv[0]*vd[1]*vv[2] , vv[0]*vv[1]*vd[2] );
#endif // GRADIENT_DOMAIN_SOLUTION
}

template< class Real >
int Octree< Real >::GetMatrixRowSize( const typename TreeOctNode::Neighbors5& neighbors5 , bool symmetric ) const
{
	int count = 0;
	int nodeIndex = neighbors5.neighbors[2][2][2]->nodeData.nodeIndex;
	const TreeOctNode* const * _nodes = &neighbors5.neighbors[0][0][0];
	if( symmetric )
	{
		for( int i=0 ; i<125 ; i++ ) if( _nodes[i] && _nodes[i]->nodeData.nodeIndex>=nodeIndex ) count++;
	}
	else
	{
		for( int i=0 ; i<125 ; i++ ) if( _nodes[i] ) count++;
	}
	return count;
}

template< class Real >
int Octree< Real >::SetMatrixRow( const PointInfo& pointInfo , const typename TreeOctNode::Neighbors5& neighbors5 , Pointer( MatrixEntry< Real > ) row , int offset , const typename BSplineData< 2 >::Integrator& integrator , const Stencil< double , 5 >& stencil , bool symmetric ) const
{
	const std::vector< _PointData >& points = pointInfo.points;
	bool hasYZPoints[3] , hasZPoints[3][3];
	Real diagonal = 0;
	Real splineValues[3*3*3*3*3];
	memset( splineValues , 0 , sizeof( Real ) * 3 * 3 * 3 * 3 * 3 );

	int count = 0;
	const TreeOctNode* node = neighbors5.neighbors[2][2][2];

	bool isInterior;
	int d , off[3];
	node->depthAndOffset( d , off );

	int o = _boundaryType==0 ? ( 1<<(d-2) ) : 0;
	int mn = 2+o , mx = (1<<d)-2-o;
	isInterior = ( off[0]>=mn && off[0]<mx && off[1]>=mn && off[1]<mx && off[2]>=mn && off[2]<mx );


	if( _constrainValues )
	{
		int idx[3] ; node->centerIndex( idx );
		for( int j=0 ; j<3 ; j++ )
		{
			hasYZPoints[j] = false;
			for( int k=0 ; k<3 ; k++ )
			{
				hasZPoints[j][k] = false;
				for( int l=0 ; l<3 ; l++ )
				{
					//neighbors5的中心三个
					const TreeOctNode* _node = neighbors5.neighbors[j+1][k+1][l+1];
					if( _node && pointInfo.pointIndex( _node )!=-1 )
					{
						const _PointData& pData = points[ pointInfo.pointIndex( _node ) ];
						//j,k,l各3维，元素为9维的
						Real* _splineValues = splineValues + 3*3*(3*(3*j+k)+l);
						Real weight = pData.weight;
						Point3D< Real > p = pData.position;
						for( int s=0 ; s<3 ; s++ )
						{
#if ROBERTO_TOLDO_FIX
							if( idx[0]+j-s>=0 && idx[0]+j-s<((2<<d)-1) ) _splineValues[3*0+s] = Real( _fData.baseBSplines[ idx[0]+j-s][s]( p[0] ) );
							if( idx[1]+k-s>=0 && idx[1]+k-s<((2<<d)-1) ) _splineValues[3*1+s] = Real( _fData.baseBSplines[ idx[1]+k-s][s]( p[1] ) );
							if( idx[2]+l-s>=0 && idx[2]+l-s<((2<<d)-1) ) _splineValues[3*2+s] = Real( _fData.baseBSplines[ idx[2]+l-s][s]( p[2] ) );
#else // !ROBERTO_TOLDO_FIX
							_splineValues[3*0+s] = Real( _fData.baseBSplines[ idx[0]+j-s][s]( p[0] ) );
							_splineValues[3*1+s] = Real( _fData.baseBSplines[ idx[1]+k-s][s]( p[1] ) );
							_splineValues[3*2+s] = Real( _fData.baseBSplines[ idx[2]+l-s][s]( p[2] ) );
#endif // ROBERTO_TOLDO_FIX
						}
						Real value = _splineValues[3*0+j] * _splineValues[3*1+k] * _splineValues[3*2+l];
						Real weightedValue = value * weight;
						for( int s=0 ; s<3 ; s++ ) _splineValues[3*0+s] *= weightedValue;
						diagonal += value * value * weight;
						hasYZPoints[j] = hasZPoints[j][k] = true;
					}
				}
			}
		}
	}

	Real pointValues[5][5][5];
	if( _constrainValues )
	{
		memset( pointValues , 0 , sizeof(Real)*5*5*5 );
		for( int i=0 ; i<3 ; i++ ) if( hasYZPoints[i] )
			for( int j=0 ; j<3 ; j++ ) if( hasZPoints[i][j] )
				for( int k=0 ; k<3 ; k++ )
				{
					const TreeOctNode* _node = neighbors5.neighbors[i+1][j+1][k+1];
					if( _node && pointInfo.pointIndex( _node )!=-1 )
						{
							const Real* _splineValuesX = splineValues + 3*(3*(3*(3*i+j)+k)+0)+2;
							const Real* _splineValuesY = splineValues + 3*(3*(3*(3*i+j)+k)+1)+2;
							const Real* _splineValuesZ = splineValues + 3*(3*(3*(3*i+j)+k)+2)+2;
							for( int ii=0 ; ii<=2 ; ii++ )
							{
								Real splineValue = _splineValuesX[-ii];
								for( int jj=0 ; jj<=2 ; jj++ )
								{
									Real* _pointValues = pointValues[i+ii][j+jj]+k;
									Real _splineValue = splineValue * _splineValuesY[-jj];
									for( int kk=0 ; kk<=2 ; kk++ ) _pointValues[kk] += _splineValue * _splineValuesZ[-kk];
								}
							}
						}
				}
	}

	pointValues[2][2][2] = diagonal;
	int nodeIndex = neighbors5.neighbors[2][2][2]->nodeData.nodeIndex;
	if( isInterior ) // General case, so try to make fast
	{
		const TreeOctNode* const * _nodes = &neighbors5.neighbors[0][0][0];
		const double* _stencil = &stencil.values[0][0][0];
		Real* _values = &pointValues[0][0][0];
		if( _constrainValues ) for( int i=0 ; i<125 ; i++ ) _values[i] = Real( _stencil[i] + _values[i] );
		else                   for( int i=0 ; i<125 ; i++ ) _values[i] = Real( _stencil[i] );
		if( symmetric ) pointValues[2][2][2] /= 2;
		//_values[5*5*2+5*2+2]即values[2][2][2]
		row[count++] = MatrixEntry< Real >( nodeIndex-offset , _values[5*5*2+5*2+2] );
		if( symmetric )
		{
			for( int i=0 ; i<125 ; i++ ) if( i!=(5*5*2+5*2+2) && _nodes[i] && _nodes[i]->nodeData.nodeIndex>=nodeIndex )
				row[count++] = MatrixEntry< Real >( _nodes[i]->nodeData.nodeIndex-offset , _values[i] );
		}
		else
		{
			for( int i=0 ; i<125 ; i++ ) if( i!=(5*5*2+5*2+2) && _nodes[i] )
				row[count++] = MatrixEntry< Real >( _nodes[i]->nodeData.nodeIndex-offset , _values[i] );
		}
	}
	else
	{
		int d , off[3];
		node->depthAndOffset( d , off );
		Real temp = Real( GetLaplacian( integrator , d , off , off , false ) );
		if( _constrainValues ) temp += pointValues[2][2][2];
		if( symmetric ) temp /= 2;
		row[count++] = MatrixEntry< Real >( nodeIndex-offset , temp );
		for( int x=0 ; x<5 ; x++ ) for( int y=0 ; y<5 ; y++ ) for( int z=0 ; z<5 ; z++ )
			if( (x!=2 || y!=2 || z!=2) && neighbors5.neighbors[x][y][z] && neighbors5.neighbors[x][y][z]->nodeData.nodeIndex>=0 && ( !symmetric || neighbors5.neighbors[x][y][z]->nodeData.nodeIndex>=nodeIndex ) )
			{
				const TreeOctNode* _node = neighbors5.neighbors[x][y][z];
				int _d , _off[3];
				_node->depthAndOffset( _d , _off );
				Real temp = Real( GetLaplacian( integrator , d , off , _off , false ) );
				if( _constrainValues ) temp += pointValues[x][y][z];
				if( symmetric && x==2 && y==2 && z==2 ) temp /= 2;
				row[count++] = MatrixEntry< Real >( _node->nodeData.nodeIndex-offset , temp );
			}
	}
	return count;
}
// if( scatter ) normals come from the center node
// else          normals come from the neighbors
template< class Real >
void Octree< Real >::SetDivergenceStencil( int depth , const typename BSplineData< 2 >::Integrator& integrator , Stencil< Point3D< double > , 5 >& stencil , bool scatter ) const
{
	if( depth<2 ) return;
	int center = 1<<(depth-1);
	int offset[] = { center , center , center };
	for( int x=0 ; x<5 ; x++ ) for( int y=0 ; y<5 ; y++ ) for( int z=0 ; z<5 ; z++ )
	{
		int _offset[] = { x+center-2 , y+center-2 , z+center-2 };
		if( scatter ) stencil.values[x][y][z] = GetDivergence1( integrator , depth , offset , _offset , false );
		else          stencil.values[x][y][z] = GetDivergence2( integrator , depth , offset , _offset , false );
	}
}
template< class Real >
void Octree< Real >::SetDivergenceStencils( int depth , const typename BSplineData< 2 >::Integrator& integrator , Stencil< Point3D< double > ,  5 > stencils[2][2][2] , bool scatter ) const
{
	if( depth<2 ) return;
	int center = 1<<(depth-1);
	for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ )
	{
		int offset[] = { center+i , center+j , center+k };
		for( int x=0 ; x<5 ; x++ ) for( int y=0 ; y<5 ; y++ ) for( int z=0 ; z<5 ; z++ )
		{
			int _offset[] = { x-2+center/2 , y-2+center/2 , z-2+center/2 };
			if( scatter ) stencils[i][j][k].values[x][y][z] = GetDivergence1( integrator , depth , offset , _offset , true );
			else          stencils[i][j][k].values[x][y][z] = GetDivergence2( integrator , depth , offset , _offset , true );
		}
	}
}
template< class Real >
void Octree< Real >::SetLaplacianStencil( int depth , const typename BSplineData< 2 >::Integrator& integrator , Stencil< double , 5 >& stencil ) const
{
	if( depth<2 ) return;
	int center = 1<<(depth-1);
	int offset[] = { center , center , center };
	for( int x=-2 ; x<=2 ; x++ ) for( int y=-2 ; y<=2 ; y++ ) for( int z=-2 ; z<=2 ; z++ )
	{
		int _offset[] = { x+center , y+center , z+center };
		stencil.values[x+2][y+2][z+2] = GetLaplacian( integrator , depth , offset , _offset , false );
	}
}
template< class Real >
void Octree< Real >::SetLaplacianStencils( int depth , const typename BSplineData< 2 >::Integrator& integrator , Stencil< double , 5 > stencils[2][2][2] ) const
{
	if( depth<2 ) return;
	int center = 1<<(depth-1);
	for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ )
	{
		int offset[] = { center+i , center+j , center+k };
		for( int x=-2 ; x<=2 ; x++ ) for( int y=-2 ; y<=2 ; y++ ) for( int z=-2 ; z<=2 ; z++ )
		{
			int _offset[] = { x+center/2 , y+center/2 , z+center/2 };
			stencils[i][j][k].values[x+2][y+2][z+2] = GetLaplacian( integrator , depth , offset , _offset , true );
		}
	}
}
template< class Real >
void Octree< Real >::SetCenterEvaluationStencil( const typename BSplineData< 2 >::template CenterEvaluator< 1 >& evaluator , int depth , Stencil< double , 3 >& stencil ) const
{
	if( depth<2 ) return;
	int center = 1<<(depth-1);
	for( int x=0 ; x<3 ; x++ ) for( int y=0 ; y<3 ; y++ ) for( int z=0 ; z<3 ; z++ )
	{
		int off[] = { center+x-1 , center+y-1 , center+z-1 };
		stencil.values[x][y][z] = Real( evaluator.value( depth , center , off[0] , false , false ) * evaluator.value( depth , center , off[1] , false , false ) * evaluator.value( depth , center , off[2] , false , false ) );
	}
}
template< class Real >
void Octree< Real >::SetCenterEvaluationStencils( const typename BSplineData< 2 >::template CenterEvaluator< 1 >& evaluator , int depth , Stencil< double , 3 > stencils[8] ) const
{
	if( depth<3 ) return;
	int center = 1<<(depth-1);
	for( int cx=0 ; cx<2 ; cx++ ) for( int cy=0 ; cy<2 ; cy++ ) for( int cz=0 ; cz<2 ; cz++ )
	{
		int idx[] = { center+cx , center+cy , center+cz };
		for( int x=0 ; x<3 ; x++ ) for( int y=0 ; y<3 ; y++ ) for( int z=0 ; z<3 ; z++ )
		{
			int off[] = { center/2+x-1 , center/2+y-1 , center/2+z-1 };
			stencils[Cube::CornerIndex( cx , cy , cz ) ].values[x][y][z] = Real( evaluator.value( depth , idx[0] , off[0] , false , true ) * evaluator.value( depth , idx[1] , off[1] , false , true ) * evaluator.value( depth , idx[2] , off[2] , false , true ) );
		}
	}
}
template< class Real >
void Octree< Real >::SetCornerEvaluationStencil( const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< double , 3 > stencil[8] ) const
{
	if( depth<2 ) return;
	int center = 1<<(depth-1);
	for( int cx=0 ; cx<2 ; cx++ ) for( int cy=0 ; cy<2 ; cy++ ) for( int cz=0 ; cz<2 ; cz++ )
	{
		int c = Cube::CornerIndex( cx , cy , cz );
		for( int x=0 ; x<3 ; x++ ) for( int y=0 ; y<3 ; y++ ) for( int z=0 ; z<3 ; z++ )
		{
			int off[] = { center+x-1 , center+y-1 , center+z-1 };
			stencil[c].values[x][y][z] = evaluator.value( depth , center , cx , off[0] , false , false ) * evaluator.value( depth , center , cy , off[1] , false , false ) * evaluator.value( depth , center , cz , off[2] , false , false );
		}
	}
}
template< class Real >
void Octree< Real >::SetCornerEvaluationStencils( const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< double , 3 > stencils[8][8] ) const
{
	if( depth<3 ) return;
	int center = 1<<(depth-1);
	for( int cx=0 ; cx<2 ; cx++ ) for( int cy=0 ; cy<2 ; cy++ ) for( int cz=0 ; cz<2 ; cz++ )
	{
		int c = Cube::CornerIndex( cx , cy , cz );
		for( int _cx=0 ; _cx<2 ; _cx++ ) for( int _cy=0 ; _cy<2 ; _cy++ ) for( int _cz=0 ; _cz<2 ; _cz++ )
		{
			int _c = Cube::CornerIndex( _cx , _cy , _cz );
			int idx[] = { center+_cx , center+_cy , center+_cz };
			for( int x=0 ; x<3 ; x++ ) for( int y=0 ; y<3 ; y++ ) for( int z=0 ; z<3 ; z++ )
			{
				int off[] = { center/2+x-1 , center/2+y-1 , center/2+z-1 };
				stencils[c][_c].values[x][y][z] = evaluator.value( depth , idx[0] , cx , off[0] , false , true ) * evaluator.value( depth , idx[1] , cy , off[1] , false , true ) * evaluator.value( depth , idx[2] , cz , off[2] , false , true );
			}
		}
	}
}
template< class Real >
void Octree< Real >::SetCornerNormalEvaluationStencil( const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< Point3D< double > , 3 > stencil[8] ) const
{
	if( depth<2 ) return;
	int center = 1<<(depth-1);
	for( int cx=0 ; cx<2 ; cx++ ) for( int cy=0 ; cy<2 ; cy++ ) for( int cz=0 ; cz<2 ; cz++ )
	{
		int c = Cube::CornerIndex( cx , cy , cz );
		for( int x=0 ; x<3 ; x++ ) for( int y=0 ; y<3 ; y++ ) for( int z=0 ; z<3 ; z++ )
		{
			int off[] = { center+x-1 , center+y-1 , center+z-1 };
			double v [] = { evaluator.value( depth , center , cx , off[0] , false , false ) , evaluator.value( depth , center , cy , off[1] , false , false ) , evaluator.value( depth , center , cz , off[2] , false , false ) };
			double dv[] = { evaluator.value( depth , center , cx , off[0] , true  , false ) , evaluator.value( depth , center , cy , off[1] , true  , false ) , evaluator.value( depth , center , cz , off[2] , true  , false ) };
			stencil[c].values[x][y][z] = Point3D< double >( dv[0]*v[1]*v[2] , v[0]*dv[1]*v[2] , v[0]*v[1]*dv[2] );
		}
	}
}
template< class Real >
void Octree< Real >::SetCornerNormalEvaluationStencils( const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< Point3D< double > , 3 > stencils[8][8] ) const
{
	if( depth<3 ) return;
	int center = 1<<(depth-1);
	for( int cx=0 ; cx<2 ; cx++ ) for( int cy=0 ; cy<2 ; cy++ ) for( int cz=0 ; cz<2 ; cz++ )
	{
		int c = Cube::CornerIndex( cx , cy , cz );	// Which corner of the finer cube
		for( int _cx=0 ; _cx<2 ; _cx++ ) for( int _cy=0 ; _cy<2 ; _cy++ ) for( int _cz=0 ; _cz<2 ; _cz++ )
		{
			int _c = Cube::CornerIndex( _cx , _cy , _cz );	// Which child node
			int idx[] = { center+_cx , center+_cy , center+_cz };
			for( int x=0 ; x<3 ; x++ ) for( int y=0 ; y<3 ; y++ ) for( int z=0 ; z<3 ; z++ )
			{
				int off[] = { center/2+x-1 , center/2+y-1 , center/2+z-1 };
				double v [] = { evaluator.value( depth , idx[0] , cx , off[0] , false , true ) , evaluator.value( depth , idx[1] , cy , off[1] , false , true ) , evaluator.value( depth , idx[2] , cz , off[2] , false , true ) };
				double dv[] = { evaluator.value( depth , idx[0] , cx , off[0] , true  , true ) , evaluator.value( depth , idx[1] , cy , off[1] , true  , true ) , evaluator.value( depth , idx[2] , cz , off[2] , true  , true ) };
				stencils[c][_c].values[x][y][z] = Point3D< double >( dv[0]*v[1]*v[2] , v[0]*dv[1]*v[2] , v[0]*v[1]*dv[2] );
			}
		}
	}
}
template< class Real >
void Octree< Real >::SetCornerNormalEvaluationStencil( const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< Point3D< double > , 5 > stencil[8] ) const
{
	if( depth<2 ) return;
	int center = 1<<(depth-1);
	for( int cx=0 ; cx<2 ; cx++ ) for( int cy=0 ; cy<2 ; cy++ ) for( int cz=0 ; cz<2 ; cz++ )
	{
		int c = Cube::CornerIndex( cx , cy , cz );
		for( int x=0 ; x<5 ; x++ ) for( int y=0 ; y<5 ; y++ ) for( int z=0 ; z<5 ; z++ )
		{
			int off[] = { center+x-2 , center+y-2 , center+z-2 };
			double v [] = { evaluator.value( depth , center , cx , off[0] , false , false ) , evaluator.value( depth , center , cy , off[1] , false , false ) , evaluator.value( depth , center , cz , off[2] , false , false ) };
			double dv[] = { evaluator.value( depth , center , cx , off[0] , true  , false ) , evaluator.value( depth , center , cy , off[1] , true  , false ) , evaluator.value( depth , center , cz , off[2] , true  , false ) };
			stencil[c].values[x][y][z] = Point3D< double >( dv[0]*v[1]*v[2] , v[0]*dv[1]*v[2] , v[0]*v[1]*dv[2] );
		}
	}
}
template< class Real >
void Octree< Real >::SetCornerNormalEvaluationStencils( const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator , int depth , Stencil< Point3D< double > , 5 > stencils[8][8] ) const
{
	if( depth<3 ) return;
	int center = 1<<(depth-1);
	for( int cx=0 ; cx<2 ; cx++ ) for( int cy=0 ; cy<2 ; cy++ ) for( int cz=0 ; cz<2 ; cz++ )
	{
		int c = Cube::CornerIndex( cx , cy , cz );	// Which corner of the finer cube
		for( int _cx=0 ; _cx<2 ; _cx++ ) for( int _cy=0 ; _cy<2 ; _cy++ ) for( int _cz=0 ; _cz<2 ; _cz++ )
		{
			int _c = Cube::CornerIndex( _cx , _cy , _cz );	// Which child node
			int idx[] = { center+_cx , center+_cy , center+_cz };
			for( int x=0 ; x<5 ; x++ ) for( int y=0 ; y<5 ; y++ ) for( int z=0 ; z<5 ; z++ )
			{
				int off[] = { center/2+x-2 , center/2+y-2 , center/2+z-2 };
				double v [] = { evaluator.value( depth , idx[0] , cx , off[0] , false , true ) , evaluator.value( depth , idx[1] , cy , off[1] , false , true ) , evaluator.value( depth , idx[2] , cz , off[2] , false , true ) };
				double dv[] = { evaluator.value( depth , idx[0] , cx , off[0] , true  , true ) , evaluator.value( depth , idx[1] , cy , off[1] , true  , true ) , evaluator.value( depth , idx[2] , cz , off[2] , true  , true ) };
				stencils[c][_c].values[x][y][z] = Point3D< double >( dv[0]*v[1]*v[2] , v[0]*dv[1]*v[2] , v[0]*v[1]*dv[2] );
			}
		}
	}
}

template< class Real >
void Octree< Real >::UpdateCoarserSupportBounds( const TreeOctNode* node , int& startX , int& endX , int& startY , int& endY , int& startZ , int& endZ )
{
	if( node->parent )
	{
		int x , y , z , c = int( node - node->parent->children );
		Cube::FactorCornerIndex( c , x , y , z );
		if( x==0 ) endX = 4;
		else     startX = 1;
		if( y==0 ) endY = 4;
		else     startY = 1;
		if( z==0 ) endZ = 4;
		else     startZ = 1;
	}
}
// Given the solution @( depth ) add to the met constraints @( depth-1 )
template< class Real >
void Octree< Real >::UpdateConstraintsFromFiner( const typename BSplineData< 2 >::Integrator& integrator , int depth , const SortedTreeNodes& sNodes , ConstPointer( Real ) fineSolution , Pointer( Real ) coarseConstraints ) const
{
	if( depth<=_minDepth ) return;
	Stencil< double , 5 > stencils[2][2][2];
	// Get the stencil describing the Laplacian relating coefficients @(depth) with coefficients @(depth-1)
	SetLaplacianStencils( depth , integrator , stencils );
	size_t start = sNodes.nodeCount[depth] , end = sNodes.nodeCount[depth+1] , range = end-start;
	int lStart = sNodes.nodeCount[depth-1];
	memset( coarseConstraints , 0 , sizeof(Real)*(sNodes.nodeCount[depth]-sNodes.nodeCount[depth-1]) );

	// Iterate over the nodes @( depth )
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys( std::max< int >( 1 , threads ) );
	for( int i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( depth-1 );
#pragma omp parallel for num_threads( threads )
	for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ )
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[ omp_get_thread_num() ];
		TreeOctNode* node = sNodes.treeNodes[i];

		bool insetSupported = _boundaryType!=0 || _IsInsetSupported( node );

		// Offset the coarser constraints using the solution from the current resolutions.
		int x , y , z , c;
		c = int( node - node->parent->children );
		Cube::FactorCornerIndex( c , x , y , z );
		if( insetSupported )
		{
			typename TreeOctNode::Neighbors5 pNeighbors5;
			neighborKey.getNeighbors( node->parent , pNeighbors5 );
			const Stencil< double , 5 >& lapStencil = stencils[x][y][z];

			Pointer( Real ) __coarseConstraints = coarseConstraints-lStart;
			bool isInterior;
			int d , off[3];
			{
				node->depthAndOffset( d , off );
				int o = _boundaryType==0 ? (1<<(d-2) ) : 0;
				int mn = 4+o , mx = (1<<d)-4-o;
				isInterior = ( off[0]>=mn && off[0]<mx && off[1]>=mn && off[1]<mx && off[2]>=mn && off[2]<mx );
			}
			// Offset the constraints using the solution from finer resolutions.
			int startX = 0 , endX = 5 , startY = 0 , endY = 5 , startZ = 0 , endZ = 5;
			UpdateCoarserSupportBounds( node , startX , endX , startY  , endY , startZ , endZ );

			Real solution = fineSolution[ node->nodeData.nodeIndex-start ];
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
				if( pNeighbors5.neighbors[x][y][z] && pNeighbors5.neighbors[x][y][z]->nodeData.nodeIndex>=0 )
				{
					const TreeOctNode* _node = pNeighbors5.neighbors[x][y][z];
					if( isInterior )
#pragma omp atomic
						__coarseConstraints[ _node->nodeData.nodeIndex ] += Real( lapStencil.values[x][y][z] * solution );
					else
					{
						int _d , _off[3];
						_node->depthAndOffset( _d , _off );
#pragma omp atomic
						__coarseConstraints[ _node->nodeData.nodeIndex ] += Real( GetLaplacian( integrator , d , off , _off , true ) * solution );
					}
				}
		}
	}
}

template< class Real >
void Octree< Real >::UpdateConstraintsFromCoarser( const PointInfo& pointInfo , const typename TreeOctNode::Neighbors5& neighbors5 , const typename TreeOctNode::Neighbors5& pNeighbors5 , TreeOctNode* node , Pointer( Real ) constraints , ConstPointer( Real ) metSolution , const typename BSplineData< 2 >::Integrator& integrator , const Stencil< double , 5 >& lapStencil ) const
{
	const std::vector< _PointData >& points = pointInfo.points;
	if( node->depth()<=_minDepth ) return;
	bool isInterior;
	int d , off[3];
	{
		node->depthAndOffset( d , off );
		int o = _boundaryType==0 ? (1<<(d-2) ) : 0;
		int mn = 4+o , mx = (1<<d)-4-o;
		isInterior = ( off[0]>=mn && off[0]<mx && off[1]>=mn && off[1]<mx && off[2]>=mn && off[2]<mx );
	}
	Real constraint = Real( 0 );
	// Offset the constraints using the solution from lower resolutions.
	int startX = 0 , endX = 5 , startY = 0 , endY = 5 , startZ = 0 , endZ = 5;
	UpdateCoarserSupportBounds( node , startX , endX , startY  , endY , startZ , endZ );

	for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
		if( pNeighbors5.neighbors[x][y][z] && pNeighbors5.neighbors[x][y][z]->nodeData.nodeIndex>=0 )
		{
			const TreeOctNode* _node = pNeighbors5.neighbors[x][y][z];
			Real _solution = metSolution[ _node->nodeData.nodeIndex ];
			{
				if( isInterior ) constraints[ node->nodeData.nodeIndex ] -= Real( lapStencil.values[x][y][z] * _solution );
				else
				{
					int _d , _off[3];
					_node->depthAndOffset( _d , _off );
					constraints[ node->nodeData.nodeIndex ] -= Real( GetLaplacian( integrator , d , off , _off , true ) * _solution );
				}
			}
		}
	if( _constrainValues )
	{
		double constraint = 0;
		int idx[3] ;
		node->centerIndex( idx );
		// Evaluate the current node's basis function at adjacent points
		for( int x=1 ; x<4 ; x++ ) for( int y=1 ; y<4 ; y++ ) for( int z=1 ; z<4 ; z++ )
			if( neighbors5.neighbors[x][y][z] && pointInfo.pointIndex( neighbors5.neighbors[x][y][z] )!=-1 )
			{
				const _PointData& pData = points[ pointInfo.pointIndex( neighbors5.neighbors[x][y][z] ) ];
				Real weightedPointValue = pData.weightedCoarserValue;
				Point3D< Real > p = pData.position;
				constraint += 
					_fData.baseBSplines[idx[0]][x-1]( p[0] ) *
					_fData.baseBSplines[idx[1]][y-1]( p[1] ) *
					_fData.baseBSplines[idx[2]][z-1]( p[2] ) * 
					weightedPointValue;
			}
		constraints[ node->nodeData.nodeIndex ] -= Real( constraint );
	}
}
struct UpSampleData
{
	int start;
	double v[2];
	UpSampleData( void ) { start = 0 , v[0] = v[1] = 0.; }
	UpSampleData( int s , double v1 , double v2 ) { start = s , v[0] = v1 , v[1] = v2; }
};
template< class Real >
template< class C >
void Octree< Real >::DownSample( int depth , const SortedTreeNodes& sNodes , ConstPointer( C ) fineConstraints , Pointer( C ) coarseConstraints ) const
{
	if( depth==0 ) return;
	double cornerValue;
	if     ( _boundaryType==-1 ) cornerValue = 0.50;
	else if( _boundaryType== 1 ) cornerValue = 1.00;
	else                         cornerValue = 0.75;
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys( std::max< int >( 1 , threads ) );
	for( int i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( depth );
#pragma omp parallel for num_threads( threads )
	for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ )
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[ omp_get_thread_num() ];
		int d , off[3];
		UpSampleData usData[3];
		sNodes.treeNodes[i]->depthAndOffset( d , off );
		for( int dd=0 ; dd<3 ; dd++ )
		{
			if     ( off[dd]  ==0          ) usData[dd] = UpSampleData( 1 , cornerValue , 0.00 );
			else if( off[dd]+1==(1<<depth) ) usData[dd] = UpSampleData( 0 , 0.00 , cornerValue );
			else if( off[dd]%2             ) usData[dd] = UpSampleData( 1 , 0.75 , 0.25 );
			else                             usData[dd] = UpSampleData( 0 , 0.25 , 0.75 );
		}
		typename TreeOctNode::Neighbors3& neighbors = neighborKey.getNeighbors( sNodes.treeNodes[i]->parent );
		C c = fineConstraints[ i-sNodes.nodeCount[depth] ];
		for( int ii=0 ; ii<2 ; ii++ )
		{
			int _ii = ii + usData[0].start;
			C cx = C( c*usData[0].v[ii] );
			for( int jj=0 ; jj<2 ; jj++ )
			{
				int _jj = jj + usData[1].start;
				C cxy = C( cx*usData[1].v[jj] );
				for( int kk=0 ; kk<2 ; kk++ )
				{
					int _kk = kk + usData[2].start;
					TreeOctNode* pNode = neighbors.neighbors[_ii][_jj][_kk];
					if( pNode )
#pragma omp atomic
						coarseConstraints[ pNode->nodeData.nodeIndex-sNodes.nodeCount[depth-1] ] += C( cxy*usData[2].v[kk] );
				}
			}
		}
	}
}
template< class Real >
template< class C >
void Octree< Real >::UpSample( int depth , const SortedTreeNodes& sNodes , ConstPointer( C ) coarseCoefficients , Pointer( C ) fineCoefficients ) const
{
	double cornerValue;
	if     ( _boundaryType==-1 ) cornerValue = 0.50;
	else if( _boundaryType== 1 ) cornerValue = 1.00;
	else                         cornerValue = 0.75;
	if( depth<=_minDepth ) return;

	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys( std::max< int >( 1 , threads ) );
	for( int i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( depth-1 );
#pragma omp parallel for num_threads( threads )
	for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ )
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[ omp_get_thread_num() ];
		bool isInterior = true;
		TreeOctNode* node = sNodes.treeNodes[i];
		int d , off[3];
		UpSampleData usData[3];
		node->depthAndOffset( d , off );
		for( int d=0 ; d<3 ; d++ )
		{
			if     ( off[d]  ==0          ) usData[d] = UpSampleData( 1 , cornerValue , 0.00 ) , isInterior = false;
			else if( off[d]+1==(1<<depth) ) usData[d] = UpSampleData( 0 , 0.00 , cornerValue ) , isInterior = false;
			else if( off[d]%2             ) usData[d] = UpSampleData( 1 , 0.75 , 0.25 );
			else                            usData[d] = UpSampleData( 0 , 0.25 , 0.75 );
		}
		typename TreeOctNode::Neighbors3& neighbors = neighborKey.getNeighbors( node->parent );
		for( int ii=0 ; ii<2 ; ii++ )
		{
			int _ii = ii + usData[0].start;
			double dx = usData[0].v[ii];
			for( int jj=0 ; jj<2 ; jj++ )
			{
				int _jj = jj + usData[1].start;
				double dxy = dx * usData[1].v[jj];
				for( int kk=0 ; kk<2 ; kk++ )
				{
					int _kk = kk + usData[2].start;
					TreeOctNode* node = neighbors.neighbors[_ii][_jj][_kk];
					if( node )
					{
						double dxyz = dxy * usData[2].v[kk];
						int _i = node->nodeData.nodeIndex;
						fineCoefficients[ i-sNodes.nodeCount[depth] ] += coarseCoefficients[ _i-sNodes.nodeCount[depth-1] ] * Real( dxyz );
					}
				}
			}
		}
	}
}
// At each point @( depth ), evaluate the met solution @( depth-1 )
template< class Real >
void Octree< Real >::SetPointValuesFromCoarser( PointInfo& pointInfo , int depth , const SortedTreeNodes& sNodes , ConstPointer( Real ) coarseCoefficients )
{
	std::vector< _PointData >& points = pointInfo.points;
	// For every node at the current depth
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys( std::max< int >( 1 , threads ) );
	for( int i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( depth );
#pragma omp parallel for num_threads( threads )
	for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ )
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[ omp_get_thread_num() ];
		int pIdx = pointInfo.pointIndex( sNodes.treeNodes[i] );
		if( pIdx!=-1 )
		{
			neighborKey.getNeighbors( sNodes.treeNodes[i] );
			points[ pIdx ].weightedCoarserValue = _WeightedCoarserFunctionValue( points[pIdx] , neighborKey , sNodes.treeNodes[i] , coarseCoefficients-_sNodes.nodeCount[depth-1] );
		}
	}
}
template< class Real >
Real Octree< Real >::_WeightedCoarserFunctionValue( const _PointData& pointData , const typename TreeOctNode::NeighborKey3& neighborKey , const TreeOctNode* pointNode , ConstPointer( Real ) coarseCoefficients ) const
{
	double pointValue = 0;
	int depth = pointNode->depth();
	if( _boundaryType==-1 && depth==0 ) return Real(-0.5) * pointData.weight;

	if( depth<=_minDepth ) return Real(0.);

	Real weight       = pointData.weight;
	Point3D< Real > p = pointData.position;

	// Iterate over all basis functions that overlap the point at the coarser resolutions
	{
		int d , _idx[3];
		const typename TreeOctNode::Neighbors3& neighbors = neighborKey.neighbors[depth-1];
		neighbors.neighbors[1][1][1]->depthAndOffset( d , _idx );
		_idx[0] = BinaryNode::CenterIndex( d , _idx[0]-1 );
		_idx[1] = BinaryNode::CenterIndex( d , _idx[1]-1 );
		_idx[2] = BinaryNode::CenterIndex( d , _idx[2]-1 );

		for( int j=0 ; j<3 ; j++ )
		{
#if ROBERTO_TOLDO_FIX
			double xValue = 0;
			if( _idx[0]+j>=0 && _idx[0]+j<((1<<depth)-1) ) xValue = _fData.baseBSplines[ _idx[0]+j ][2-j]( p[0] );
			else continue;
#else // !ROBERTO_TOLDO_FIX
			double xValue = _fData.baseBSplines[ _idx[0]+j ][2-j]( p[0] );
#endif // ROBERTO_TOLDO_FIX
			for( int k=0 ; k<3 ; k++ )
			{
#if ROBERTO_TOLDO_FIX
				double xyValue = 0;
				if( _idx[1]+k>=0 && _idx[1]+k<((1<<depth)-1) ) xyValue = xValue * _fData.baseBSplines[ _idx[1]+k ][2-k]( p[1] );
				else continue;
#else // !ROBERTO_TOLDO_FIX
				double xyValue = xValue * _fData.baseBSplines[ _idx[1]+k ][2-k]( p[1] );
#endif // ROBERTO_TOLDO_FIX
				double _pointValue = 0;
				for( int l=0 ; l<3 ; l++ )
				{
					const TreeOctNode* basisNode = neighbors.neighbors[j][k][l];
#if ROBERTO_TOLDO_FIX
					if( basisNode && basisNode->nodeData.nodeIndex>=0 && _idx[2]+l>=0 && _idx[2]+l<((1<<depth)-1) )
						_pointValue += _fData.baseBSplines[ _idx[2]+l ][2-l]( p[2] ) * double( coarseCoefficients[basisNode->nodeData.nodeIndex] );
#else // !ROBERTO_TOLDO_FIX
					if( basisNode && basisNode->nodeData.nodeIndex>=0 )
						_pointValue += _fData.baseBSplines[ _idx[2]+l ][2-l]( p[2] ) * double( coarseCoefficients[basisNode->nodeData.nodeIndex] );
#endif // ROBERTO_TOLDO_FIX
				}
				pointValue += _pointValue * xyValue;
			}
		}
	}
	if( _boundaryType==-1 ) pointValue -= 0.5;
	return Real( pointValue * weight );
}
template< class Real >
void Octree< Real >::SetPointConstraintsFromFiner( const PointInfo& pointInfo , int depth , const SortedTreeNodes& sNodes , ConstPointer( Real ) finerCoefficients , Pointer( Real ) coarserConstraints ) const
{
	const std::vector< _PointData >& points = pointInfo.points;
	// Note: We can't iterate over the finer point nodes as the point weights might be
	// scaled incorrectly, due to the adaptive exponent. So instead, we will iterate
	// over the coarser nodes and evaluate the finer solution at the associated points.
	if( !depth ) return;
	size_t start = sNodes.nodeCount[depth-1] , end = sNodes.nodeCount[depth] , range = end-start;
	memset( coarserConstraints , 0 , sizeof( Real ) * ( sNodes.nodeCount[depth]-sNodes.nodeCount[depth-1] ) );
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys( std::max< int >( 1 , threads ) );
	for( int i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( depth-1 );
#pragma omp parallel for num_threads( threads )
	for( int i=sNodes.nodeCount[depth-1] ; i<sNodes.nodeCount[depth] ; i++ )
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[ omp_get_thread_num() ];
		int pIdx = pointInfo.pointIndex( sNodes.treeNodes[i] );
		if( pIdx!=-1 )
		{
			typename TreeOctNode::Neighbors3& neighbors = neighborKey.getNeighbors( sNodes.treeNodes[i] );
			// Evaluate the solution @( depth ) at the current point @( depth-1 )
			{
				Real finerPointValue = _WeightedFinerFunctionValue( points[pIdx] , neighborKey , sNodes.treeNodes[i] , finerCoefficients-sNodes.nodeCount[depth] );
				Point3D< Real > p = points[ pIdx ].position;
				// Update constraints for all nodes @( depth-1 ) that overlap the point
				int d , idx[3];
				neighbors.neighbors[1][1][1]->depthAndOffset( d, idx );
				// Set the (offset) index to the top-left-front corner of the 3x3x3 block of b-splines
				// overlapping the point.
				idx[0] = BinaryNode::CenterIndex( d , idx[0]-1 );
				idx[1] = BinaryNode::CenterIndex( d , idx[1]-1 );
				idx[2] = BinaryNode::CenterIndex( d , idx[2]-1 );
				for( int x=0 ; x<3 ; x++ ) for( int y=0 ; y<3 ; y++ ) for( int z=0 ; z<3 ; z++ )
					if( neighbors.neighbors[x][y][z] )
					{
#pragma omp atomic
						coarserConstraints[ neighbors.neighbors[x][y][z]->nodeData.nodeIndex - sNodes.nodeCount[depth-1] ] +=
							Real(
							_fData.baseBSplines[idx[0]+x][2-x]( p[0] ) *
							_fData.baseBSplines[idx[1]+y][2-y]( p[1] ) *
							_fData.baseBSplines[idx[2]+z][2-z]( p[2] ) * 
							finerPointValue
							);
					}
			}
		}
	}
}
template< class Real >
Real Octree< Real >::_WeightedFinerFunctionValue( const _PointData& pointData , const typename TreeOctNode::NeighborKey3& neighborKey , const TreeOctNode* pointNode , ConstPointer( Real ) finerCoefficients ) const
{
	typename TreeOctNode::Neighbors3 childNeighbors;
	double pointValue = 0;
	int depth = pointNode->depth();
	Real weight       = pointData.weight;
	Point3D< Real > p = pointData.position;
	neighborKey.getChildNeighbors( p , depth , childNeighbors );
	// Iterate over all finer basis functions that overlap the point at the coarser resolutions
	int d , idx[3];
	{
		Point3D< Real > c;
		Real w;
		neighborKey.neighbors[depth].neighbors[1][1][1]->depthAndOffset( d , idx );
		neighborKey.neighbors[depth].neighbors[1][1][1]->centerAndWidth( c , w );
		d++;
		idx[0] *= 2 , idx[1] *= 2 , idx[2] *= 2;
		int cIndex=TreeOctNode::CornerIndex( c , p );
		if( cIndex&1 ) idx[0]++;
		if( cIndex&2 ) idx[1]++;
		if( cIndex&4 ) idx[2]++;
	}
	// Center the indexing at the top-left-front corner
	idx[0] = BinaryNode::CenterIndex( d , idx[0]-1 );
	idx[1] = BinaryNode::CenterIndex( d , idx[1]-1 );
	idx[2] = BinaryNode::CenterIndex( d , idx[2]-1 );

	for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ ) for( int l=0 ; l<3 ; l++ )
	{
		const TreeOctNode* basisNode = childNeighbors.neighbors[j][k][l];
		if( basisNode && basisNode->nodeData.nodeIndex>=0 )
			pointValue += 
			_fData.baseBSplines[ idx[0]+j ][2-j]( p[0] ) *
			_fData.baseBSplines[ idx[1]+k ][2-k]( p[1] ) *
			_fData.baseBSplines[ idx[2]+l ][2-l]( p[2] ) *
			double( finerCoefficients[ basisNode->nodeData.nodeIndex ] );
	}
	if( _boundaryType==-1 ) pointValue -= Real(0.5);
	return Real( pointValue * weight );
}
template< class Real >
int Octree< Real >::GetSliceMatrixAndUpdateConstraints( const PointInfo& pointInfo , SparseMatrix< Real >& matrix , Pointer( Real ) constraints , const typename BSplineData< 2 >::Integrator& integrator , int depth , const SortedTreeNodes& sNodes , ConstPointer( Real ) metSolution , bool coarseToFine , int nStart , int nEnd )
{
	size_t range = nEnd-nStart;
	Stencil< double , 5 > stencil , stencils[2][2][2];
	SetLaplacianStencil ( depth , integrator , stencil );
	SetLaplacianStencils( depth , integrator , stencils );
	matrix.Resize( (int)range );
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys( std::max< int >( 1 , threads ) );
	for( int i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( depth );
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<range ; i++ )
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[ omp_get_thread_num() ];
		TreeOctNode* node = sNodes.treeNodes[i+nStart];
		// Get the matrix row size
		bool insetSupported = _boundaryType!=0 || _IsInsetSupported( node );
		typename TreeOctNode::Neighbors5 neighbors5;
		if( insetSupported ) neighborKey.getNeighbors( node , neighbors5 );
		int count = insetSupported ? GetMatrixRowSize( neighbors5 , false ) : 1;

		// Allocate memory for the row
#pragma omp critical (matrix_set_row_size)
		{
			matrix.SetRowSize( i , count );
		}

		// Set the row entries
		if( insetSupported ) matrix.rowSizes[i] = SetMatrixRow( pointInfo , neighbors5 , matrix[i] , sNodes.nodeCount[depth] , integrator , stencil , false );
		else
		{
			matrix[i][0] = MatrixEntry< Real >( i , Real(1) );
			matrix.rowSizes[i] = 1;
		}

		if( depth>_minDepth )
		{
			// Offset the constraints using the solution from lower resolutions.
			int x , y , z , c;
			if( node->parent )
			{
				c = int( node - node->parent->children );
				Cube::FactorCornerIndex( c , x , y , z );
			}
			else x = y = z = 0;
			if( insetSupported && coarseToFine )
			{
				typename TreeOctNode::Neighbors5 pNeighbors5;
				neighborKey.getNeighbors( node->parent , pNeighbors5 );
				UpdateConstraintsFromCoarser( pointInfo , neighbors5 , pNeighbors5 , node , constraints , metSolution , integrator , stencils[x][y][z] );
			}
		}
	}
	return 1;
}
template< class Real >
int Octree< Real >::GetMatrixAndUpdateConstraints( const PointInfo& pointInfo , SparseSymmetricMatrix< Real >& matrix , Pointer( Real ) constraints , const typename BSplineData< 2 >::Integrator& integrator , int depth , const SortedTreeNodes& sNodes , ConstPointer( Real ) metSolution , bool coarseToFine )
{
	size_t start = sNodes.nodeCount[depth] , end = sNodes.nodeCount[depth+1] , range = end-start;
	Stencil< double , 5 > stencil , stencils[2][2][2];
	SetLaplacianStencil ( depth , integrator , stencil );
	SetLaplacianStencils( depth , integrator , stencils );
	matrix.Resize( (int)range );
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys( std::max< int >( 1 , threads ) );
	for( int i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( depth );
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<range ; i++ )
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[ omp_get_thread_num() ];
		TreeOctNode* node = sNodes.treeNodes[i+start];
		// Get the matrix row size
		bool insetSupported = _boundaryType!=0 || _IsInsetSupported( node );
		typename TreeOctNode::Neighbors5 neighbors5;
		if( insetSupported ) neighborKey.getNeighbors( node , neighbors5 );
		int count = insetSupported ? GetMatrixRowSize( neighbors5 , true ) : 1;

		// Allocate memory for the row
#pragma omp critical (matrix_set_row_size)
		matrix.SetRowSize( i , count );

		// Set the row entries
		if( insetSupported ) matrix.rowSizes[i] = SetMatrixRow( pointInfo , neighbors5 , matrix[i] , (int)start , integrator , stencil , true );
		else
		{
			matrix[i][0] = MatrixEntry< Real >( i , Real(1) );
			matrix.rowSizes[i] = 1;
		}
		if( depth>_minDepth )
		{
			// Offset the constraints using the solution from lower resolutions.
			int x , y , z , c;
			if( node->parent )
			{
				c = int( node - node->parent->children );
				Cube::FactorCornerIndex( c , x , y , z );
			}
			else x = y = z = 0;
			if( insetSupported && coarseToFine )
			{
				typename TreeOctNode::Neighbors5 pNeighbors5;
				neighborKey.getNeighbors( node->parent , pNeighbors5 );
				UpdateConstraintsFromCoarser( pointInfo , neighbors5 , pNeighbors5 , node , constraints , metSolution , integrator , stencils[x][y][z] );
			}
		}
	}
	return 1;
}

template< class Real >
Pointer( Real ) Octree< Real >::SolveSystem( PointInfo& pointInfo , Pointer( Real ) constraints , bool showResidual , int iters , int maxSolveDepth , int cgDepth , double accuracy )
{
	int iter=0;
	typename BSplineData< 2 >::Integrator integrator;
	_fData.setIntegrator( integrator , _boundaryType==0 );
	iters = std::max< int >( 0 , iters );
	if( _boundaryType==0 ) maxSolveDepth++ , cgDepth++;

	Pointer( Real ) solution = AllocPointer< Real >( _sNodes.nodeCount[_sNodes.maxDepth] );
	memset( solution , 0 , sizeof(Real)*_sNodes.nodeCount[_sNodes.maxDepth] );

	solution[0] = 0;

	//metSolution是Poisson Algorithm 2中的coarser solution
	std::vector< Real > metSolution( _sNodes.nodeCount[ _sNodes.maxDepth-1 ] , 0 );
	bool coarse_to_fine = true;
	//从min depth到max depth
	for( int d=_minDepth ; d<_sNodes.maxDepth ; d++ )
	{
		DumpOutput( "Depth[%d/%d]: %d\n" , _boundaryType==0 ? d-1 : d , _boundaryType==0 ? _sNodes.maxDepth-2 : _sNodes.maxDepth-1 , _sNodes.nodeCount[d+1]-_sNodes.nodeCount[d] );
		if( d==_minDepth )
			_SolveSystemCG( pointInfo , d , integrator , _sNodes , solution , constraints , GetPointer( metSolution ) , _sNodes.nodeCount[_minDepth+1]-_sNodes.nodeCount[_minDepth] , coarse_to_fine, showResidual, NULL , NULL , NULL );
		else
		{
			if( d>cgDepth ) iter += _SolveSystemGS( pointInfo , d , integrator , _sNodes , solution , constraints , GetPointer( metSolution ) , d>maxSolveDepth ? 0 : iters , coarse_to_fine, showResidual , NULL , NULL , NULL );
			else            iter += _SolveSystemCG( pointInfo , d , integrator , _sNodes , solution , constraints , GetPointer( metSolution ) , d>maxSolveDepth ? 0 : iters , coarse_to_fine, showResidual , NULL , NULL , NULL , accuracy );
		}
	}

	return solution;
}
template< class Real >
void Octree< Real >::_setMultiColorIndices( int start , int end , std::vector< std::vector< int > >& indices ) const
{
	const int modulus = 3;
	indices.resize( modulus*modulus*modulus );
	int count[modulus*modulus*modulus];
	memset( count , 0 , sizeof(int)*modulus*modulus*modulus );
#pragma omp parallel for num_threads( threads )
	for( int i=start ; i<end ; i++ )
	{
		int d , off[3];
		_sNodes.treeNodes[i]->depthAndOffset( d , off );
		int idx = (modulus*modulus) * ( off[2]%modulus ) + modulus * ( off[1]%modulus ) + ( off[0]%modulus );
#pragma omp atomic
		count[idx]++;
	}

	for( int i=0 ; i<modulus*modulus*modulus ; i++ ) indices[i].reserve( count[i] ) , count[i]=0;

	for( int i=start ; i<end ; i++ )
	{
		int d , off[3];
		_sNodes.treeNodes[i]->depthAndOffset( d , off );
		int idx = (modulus*modulus) * ( off[2]%modulus ) + modulus * ( off[1]%modulus ) + ( off[0]%modulus );
		indices[idx].push_back( _sNodes.treeNodes[i]->nodeData.nodeIndex - start );
	}
}
template< class Real >
int Octree< Real >::_SolveSystemGS( PointInfo& pointInfo , int depth , const typename BSplineData< 2 >::Integrator& integrator , const SortedTreeNodes& sNodes , Pointer( Real ) solution , Pointer( Real ) constraints , Pointer( Real ) metSolutionConstraints , int iters , bool coarseToFine , bool showResidual , double* bNorm2 , double* inRNorm2 , double* outRNorm2 , bool forceSilent )
{
	Pointer( Real ) metSolution = NullPointer< Real >();
	Pointer( Real ) metConstraints = NullPointer< Real >();
	if( coarseToFine ) metSolution    = metSolutionConstraints;	// This stores the up-sampled solution up to depth-2
	else               metConstraints = metSolutionConstraints; // This stores the down-sampled constraints up to depth

	double _maxMemoryUsage = maxMemoryUsage;
	maxMemoryUsage = 0;
	Vector< Real > X , B;
	int slices = 1<<depth;
	double systemTime=0. , solveTime=0. , updateTime=0. ,  evaluateTime = 0.;
	std::vector< int > offsets( slices+1 , 0 );
	for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ )
	{
		int d , off[3];
		sNodes.treeNodes[i]->depthAndOffset( d , off );
		offsets[ off[2] ]++;
	}
	for( int i=1 ; i<slices ; i++ )  offsets[i] += offsets[i-1];
	for( int i=slices ; i>=1 ; i-- ) offsets[i]  = offsets[i-1];
	offsets[0] = 0;

	X.Resize( sNodes.nodeCount[depth+1]-sNodes.nodeCount[depth] );
	B.Resize( sNodes.nodeCount[depth+1]-sNodes.nodeCount[depth] );
	if( coarseToFine )
	{
		if( depth>_minDepth )
		{
			// Up-sample the cumulative change in solution @(depth-2) into the cumulative change in solution @(depth-1)
			if( depth-2>=_minDepth ) UpSample( depth-1 , sNodes , ( ConstPointer( Real ) )metSolution+_sNodes.nodeCount[depth-2] , metSolution+_sNodes.nodeCount[depth-1] );
			// Add in the change in solution @(depth-1)
#pragma omp parallel for num_threads( threads )
			for( int i=_sNodes.nodeCount[depth-1] ; i<_sNodes.nodeCount[depth] ; i++ ) metSolution[i] += solution[i];
			// Evaluate the points @(depth) using the cumulative change in solution @(depth-1)
			if( _constrainValues )
			{
				evaluateTime = Time();
				SetPointValuesFromCoarser( pointInfo , depth , sNodes , metSolution+_sNodes.nodeCount[depth-1] );
				evaluateTime = Time() - evaluateTime;
			}
		}
	}
	else if( depth<_sNodes.maxDepth-1 )
		for( int i=_sNodes.nodeCount[depth] ; i<_sNodes.nodeCount[depth+1] ; i++ ) constraints[i] -= metConstraints[i];
	// Initialize with the previously computed solution
#pragma omp parallel for num_threads( threads )
	for( int i=_sNodes.nodeCount[depth] ; i<_sNodes.nodeCount[depth+1] ; i++ ) X[ i-_sNodes.nodeCount[depth] ] = solution[i];
	double bNorm=0 , inRNorm=0 , outRNorm=0;
	if( depth>=_minDepth )
	{
		int frontOffset = ( showResidual || inRNorm2 ) ? 2 : 0;
		int backOffset = ( showResidual || outRNorm2 ) ? 2 : 0;
		int solveSlices = std::min< int >( 2*iters-1 , slices ) , matrixSlices = std::max< int >( 1 , std::min< int >( solveSlices+frontOffset+backOffset , slices ) );
		std::vector< SparseMatrix< Real > > _M( matrixSlices );
		std::vector< std::vector< std::vector< int > > > __mcIndices( std::max< int >( 0 , solveSlices ) );

		int dir = coarseToFine ? -1 : 1 , start = coarseToFine ? slices-1 : 0 , end = coarseToFine ? -1 : slices;
		for( int frontSlice=start-frontOffset*dir , backSlice = frontSlice-2*(iters-1)*dir ; backSlice!=end+backOffset*dir ; frontSlice+=dir , backSlice+=dir )
		{
			double t;
			if( frontSlice+frontOffset*dir>=0 && frontSlice+frontOffset*dir<slices )
			{
				int s = frontSlice+frontOffset*dir , _s = s % matrixSlices;
				t = Time();
				GetSliceMatrixAndUpdateConstraints( pointInfo , _M[_s] , constraints , integrator , depth , sNodes , metSolution , coarseToFine , offsets[s]+sNodes.nodeCount[depth] , offsets[s+1]+sNodes.nodeCount[depth] );
				systemTime += Time()-t;
				Pointer( TreeOctNode* ) const nodes = sNodes.treeNodes + sNodes.nodeCount[depth];
				for( int i=offsets[s] ; i<offsets[s+1] ; i++ )
				{
					if( _boundaryType!=0 || _IsInsetSupported( nodes[i] ) ) B[i] = constraints[ nodes[i]->nodeData.nodeIndex ];
					else                                                    B[i] = Real(0);
				}
				if( showResidual || inRNorm2 )
#pragma omp parallel for num_threads( threads ) reduction( + : bNorm , inRNorm )
					for( int j=0 ; j<_M[_s].rows ; j++ )
					{
						Real temp = Real(0);
						ConstPointer( MatrixEntry< Real > ) start = _M[_s][j];
						ConstPointer( MatrixEntry< Real > ) end = start + _M[_s].rowSizes[j];
						ConstPointer( MatrixEntry< Real > ) e;
						for( e=start ; e!=end ; e++ ) temp += X[ e->N ] * e->Value;
						Real b = B[ j + offsets[s] ] ;
						bNorm += b*b;
						inRNorm += (temp-b) * (temp-b);
					}
				else if( bNorm2 )
#pragma omp parallel for num_threads( threads ) reduction( + : bNorm )
					for( int j=0 ; j<_M[_s].rows ; j++ )
					{
						Real b = B[ j + offsets[s] ] ;
						bNorm += b*b;
					}
			}
			t = Time();
			if( iters && frontSlice>=0 && frontSlice<slices )
			{
				int s = frontSlice , _s = s % matrixSlices , __s = s % solveSlices;
				for( int i=0 ; i<int( __mcIndices[__s].size() ) ; i++ ) __mcIndices[__s][i].clear();
				_setMultiColorIndices( sNodes.nodeCount[depth]+offsets[s] , sNodes.nodeCount[depth]+offsets[s+1] , __mcIndices[__s] );
			}
			for( int slice=frontSlice ; slice*dir>=backSlice*dir ; slice-=2*dir )
				if( slice>=0 && slice<slices )
				{
					int s = slice , _s = s % matrixSlices , __s = s % solveSlices;
					SparseMatrix< Real >::SolveGS( __mcIndices[__s] , _M[_s] , B , X , !coarseToFine , threads , offsets[s] );
				}
			solveTime += Time() - t;
			if( (showResidual || outRNorm2) && backSlice-backOffset*dir>=0 && backSlice-backOffset*dir<slices )
			{
				int s = backSlice-backOffset*dir , _s = s % matrixSlices;
#pragma omp parallel for num_threads( threads ) reduction( + : outRNorm )
				for( int j=0 ; j<_M[_s].rows ; j++ )
				{
					Real temp = Real(0);
					ConstPointer( MatrixEntry< Real > ) start = _M[_s][j];
					ConstPointer( MatrixEntry< Real > ) end = start + _M[_s].rowSizes[j];
					ConstPointer( MatrixEntry< Real > ) e;
					for( e=start ; e!=end ; e++ ) temp += X[ e->N ] * e->Value;
					Real b = B[ j + offsets[s] ];
					outRNorm += (temp-b) * (temp-b);
				}
			}
		}
	}

	if( bNorm2 ) bNorm2[depth] = bNorm;
	if( inRNorm2 ) inRNorm2[depth] = inRNorm;
	if( outRNorm2 ) outRNorm2[depth] = outRNorm;
	if( showResidual && iters )
	{
		for( int i=0 ; i<depth ; i++ ) printf( "  " );
		printf( "GS: %.4e -> %.4e -> %.4e (%.2e) [%d]\n" , sqrt( bNorm ) , sqrt( inRNorm ) , sqrt( outRNorm ) , sqrt( outRNorm/bNorm ) , iters );
	}

	// Copy the old solution into the buffer, write in the new solution, compute the change, and update the met constraints
#pragma omp parallel for num_threads( threads )
	for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ ) solution[i] = X[ i-sNodes.nodeCount[depth] ];
	if( !coarseToFine && depth>_minDepth )
	{
		// Explicitly compute the restriction of the met solution onto the coarser nodes
		// and down-sample the previous accumulation
		{
			UpdateConstraintsFromFiner( integrator , depth , sNodes , GetPointer( X ) , metConstraints+sNodes.nodeCount[depth-1] );
			if( _constrainValues ) SetPointConstraintsFromFiner( pointInfo , depth , sNodes , GetPointer( X ) , metConstraints+sNodes.nodeCount[depth-1] );
			if( depth<sNodes.maxDepth-1 ) DownSample( depth , sNodes , ( ConstPointer( Real ) )metConstraints+sNodes.nodeCount[depth] , metConstraints+sNodes.nodeCount[depth-1] );
		}
	}

	MemoryUsage();
	if( !forceSilent ) DumpOutput( "\tEvaluated / Got / Solved in: %6.3f / %6.3f / %6.3f\t(%.3f MB)\n" , evaluateTime , systemTime , solveTime , float( maxMemoryUsage ) );
	maxMemoryUsage = std::max< double >( maxMemoryUsage , _maxMemoryUsage );

	return iters;
}
template< class Real >
int Octree< Real >::_SolveSystemCG( PointInfo& pointInfo , int depth , const typename BSplineData< 2 >::Integrator& integrator , const SortedTreeNodes& sNodes , Pointer( Real ) solution , Pointer( Real ) constraints , Pointer( Real ) metSolutionConstraints , int iters , bool coarseToFine , bool showResidual , double* bNorm2 , double* inRNorm2 , double* outRNorm2 , double accuracy )
{
	Pointer( Real ) metSolution = NullPointer< Real >();
	Pointer( Real ) metConstraints = NullPointer< Real >();
	if( coarseToFine ) metSolution    = metSolutionConstraints;	// This stores the up-sampled solution up to depth-2
	else               metConstraints = metSolutionConstraints; // This stores the down-sampled constraints up to depth
	//metSolution存放up-sampled solution up to depth-2
	double _maxMemoryUsage = maxMemoryUsage;
	maxMemoryUsage = 0;
	int iter = 0;
	Vector< Real > X , B;
	SparseSymmetricMatrix< Real > M;
	double systemTime=0. , solveTime=0. , updateTime=0. ,  evaluateTime = 0.;
	X.Resize( sNodes.nodeCount[depth+1]-sNodes.nodeCount[depth] );
	if( coarseToFine )
	{
		if( depth>_minDepth )
		{
			// Up-sample the cumulative change in solution @(depth-2) into the cumulative change in solution @(depth-1)
			if( depth-2>=_minDepth ) UpSample( depth-1 , sNodes , ( ConstPointer( Real ) )metSolution+_sNodes.nodeCount[depth-2] , metSolution+_sNodes.nodeCount[depth-1] );
			// Add in the change in solution @(depth-1)
#pragma omp parallel for num_threads( threads )
			for( int i=_sNodes.nodeCount[depth-1] ; i<_sNodes.nodeCount[depth] ; i++ ) metSolution[i] += solution[i];
			//_constrainValues != 0说明Poisson中的point constrains系数不为0
			// Evaluate the points @(depth) using the cumulative change in solution @(depth-1),
			if( _constrainValues )
			{
				evaluateTime = Time();
				//用于求Point data中的weightedCoarserValue
				SetPointValuesFromCoarser( pointInfo , depth , sNodes , metSolution+_sNodes.nodeCount[depth-1] );
				evaluateTime = Time() - evaluateTime;
			}
		}
	}
	else if( depth<_sNodes.maxDepth-1 )
		for( int i=_sNodes.nodeCount[depth] ; i<_sNodes.nodeCount[depth+1] ; i++ ) constraints[i] -= metConstraints[i];
	// Initialize with the previously computed solution
#pragma omp parallel for num_threads( threads )
	for( int i=_sNodes.nodeCount[depth] ; i<_sNodes.nodeCount[depth+1] ; i++ ) X[ i-_sNodes.nodeCount[depth] ] = solution[i];
	systemTime = Time();
	{
		// Get the system matrix (and adjust the right-hand-side based on the coarser solution if prolonging)
		if( coarseToFine ) GetMatrixAndUpdateConstraints( pointInfo , M , constraints , integrator , depth , sNodes , metSolution           , true  );
		else               GetMatrixAndUpdateConstraints( pointInfo , M , constraints , integrator , depth , sNodes , NullPointer< Real >() , false );
		// Set the constraint vector
		B.Resize( sNodes.nodeCount[depth+1]-sNodes.nodeCount[depth] );
		for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ )
			if( _boundaryType!=0 || _IsInsetSupported( sNodes.treeNodes[i] ) ) B[i-sNodes.nodeCount[depth]] = constraints[i];
			else                                                               B[i-sNodes.nodeCount[depth]] = Real(0);
	}
	systemTime = Time()-systemTime;

	solveTime = Time();
	// Solve the linear system
	accuracy = Real( accuracy / 100000 ) * M.rows;
	int res = 1<<depth;

	MapReduceVector< Real > mrVector;
	mrVector.resize( threads , M.rows );
	bool addDCTerm = (M.rows==res*res*res && !_constrainValues && _boundaryType!=-1);
	double bNorm , inRNorm , outRNorm;
	if( showResidual || bNorm2 ) bNorm = B.Norm( 2 );
	if( showResidual || inRNorm2 ) inRNorm = ( addDCTerm ? ( B - M * X - X.Average() ) : ( B - M * X ) ).Norm( 2 );

	if( _boundaryType==0 && depth>3 ) res -= 1<<(depth-2);
	if( iters ) iter += SparseSymmetricMatrix< Real >::SolveCG( M , B , iters , X , mrVector , Real( accuracy ) , 0 , addDCTerm );
	solveTime = Time()-solveTime;
	if( showResidual || outRNorm2 ) outRNorm = ( addDCTerm ? ( B - M * X - X.Average() ) : ( B - M * X ) ).Norm( 2 );
	if( bNorm2 ) bNorm2[depth] = bNorm * bNorm;
	if( inRNorm2 ) inRNorm2[depth] = inRNorm * inRNorm;
	if( outRNorm2 ) outRNorm2[depth] = outRNorm * outRNorm;
	if( showResidual && iters )
	{
		for( int i=0 ; i<depth ; i++ ) printf( "  " );
		printf( "CG: %.4e -> %.4e -> %.4e (%.2e) [%d]\n" , bNorm , inRNorm , outRNorm , outRNorm/bNorm , iter );
	}

	// Copy the old solution into the buffer, write in the new solution, compute the change, and update the met solution
	{
#pragma omp parallel for num_threads( threads )
		for( int i=sNodes.nodeCount[depth] ; i<sNodes.nodeCount[depth+1] ; i++ ) solution[i] = X[ i-sNodes.nodeCount[depth] ];
		if( !coarseToFine && depth>_minDepth )
		{
			// Explicitly compute the restriction of the met solution onto the coarser nodes
			// and down-sample the previous accumulation
			{
				UpdateConstraintsFromFiner( integrator , depth , sNodes , GetPointer( X ) , metConstraints + sNodes.nodeCount[depth-1] );
				if( _constrainValues ) SetPointConstraintsFromFiner( pointInfo , depth , sNodes , GetPointer( X ) , metConstraints+sNodes.nodeCount[depth-1] );
				if( depth<sNodes.maxDepth-1 ) DownSample( depth , sNodes , ( ConstPointer( Real ) )metConstraints+sNodes.nodeCount[depth] , metConstraints+sNodes.nodeCount[depth-1] );
			}
		}
	}

	MemoryUsage();
	DumpOutput( "\tEvaluated / Got / Solved in: %6.3f / %6.3f / %6.3f\t(%.3f MB)\n" , evaluateTime , systemTime , solveTime , float( maxMemoryUsage ) );
	maxMemoryUsage = std::max< double >( maxMemoryUsage , _maxMemoryUsage );
	return iter;
}
template< class Real >
int Octree< Real >::HasNormals( TreeOctNode* node , const NormalInfo& normalInfo )
{
	int idx = normalInfo.normalIndex( node );
	if( idx>=0 )
	{
		const Point3D< Real >& normal = normalInfo.normals[ idx ];
		if( normal[0]!=0 || normal[1]!=0 || normal[2]!=0 ) return 1;
	}
	if( node->children ) for( int i=0 ; i<Cube::CORNERS ; i++ ) if( HasNormals( &node->children[i] , normalInfo ) ) return 1;
	return 0;
}
template< class Real >
Pointer( Real ) Octree< Real >::SetLaplacianConstraints( const NormalInfo& normalInfo )
{
	// To set the Laplacian constraints, we iterate over the
	// splatted normals and compute the dot-product of the
	// divergence of the normal field with all the basis functions.
	// Within the same depth: set directly as a gather
	// Coarser depths 
	typename BSplineData< 2 >::Integrator integrator;
	_fData.setIntegrator( integrator , _boundaryType==0 );
	int maxDepth = _sNodes.maxDepth-1;
	Point3D< Real > zeroPoint;
	zeroPoint[0] = zeroPoint[1] = zeroPoint[2] = 0;
	Pointer( Real ) constraints = AllocPointer< Real >( _sNodes.nodeCount[_sNodes.maxDepth] );
	if( !constraints ) fprintf( stderr , "[ERROR] Failed to allocate constraints: %d * %zu\n" , _sNodes.nodeCount[_sNodes.maxDepth] , sizeof( Real ) ) , exit( 0 );
	memset( constraints , 0 , sizeof(Real)*_sNodes.nodeCount[_sNodes.maxDepth] );
	Pointer( Real ) _constraints = AllocPointer< Real >( _sNodes.nodeCount[maxDepth] );
	if( !_constraints ) fprintf( stderr , "[ERROR] Failed to allocate _constraints: %d * %zu\n" , _sNodes.nodeCount[maxDepth] , sizeof( Real ) ) , exit( 0 );
	memset( _constraints , 0 , sizeof(Real)*_sNodes.nodeCount[maxDepth] );
	MemoryUsage();

	for( int d=maxDepth ; d>=(_boundaryType==0?2:0) ; d-- )
	{
		int offset = d>0 ? _sNodes.treeNodes[ _sNodes.nodeCount[d-1] ]->nodeData.nodeIndex : 0;
		Stencil< Point3D< double > , 5 > stencil , stencils[2][2][2];
		//stencil的值是某一深度两个相邻结点B样条函数之间grad(Bi)和Bj的各维度内积。stencils是相差一层的父子结点之间积分。
		SetDivergenceStencil ( d , integrator , stencil , false );
		SetDivergenceStencils( d , integrator , stencils , true );

		std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys( std::max< int >( 1 , threads ) );
		for( int i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _fData.depth );
#pragma omp parallel for num_threads( threads )
		for( int i=_sNodes.nodeCount[d] ; i<_sNodes.nodeCount[d+1] ; i++ )
		{
			typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[ omp_get_thread_num() ];
			TreeOctNode* node = _sNodes.treeNodes[i];
			int startX=0 , endX=5 , startY=0 , endY=5 , startZ=0 , endZ=5;
			int depth = node->depth();
			typename TreeOctNode::Neighbors5 neighbors5;
			neighborKey.getNeighbors( node , neighbors5 );

			bool isInterior , isInterior2;
			{
				int d , off[3];
				node->depthAndOffset( d , off );
				int o = _boundaryType==0 ? (1<<(d-2)) : 0;
				int mn = 2+o , mx = (1<<d)-2-o;
				isInterior  = ( off[0]>=mn && off[0]<mx && off[1]>=mn && off[1]<mx && off[2]>=mn && off[2]<mx );
				mn += 2 , mx -= 2;
				isInterior2 = ( off[0]>=mn && off[0]<mx && off[1]>=mn && off[1]<mx && off[2]>=mn && off[2]<mx );
			}
			int cx , cy , cz;
			if( d )
			{
				int c = int( node - node->parent->children );
				Cube::FactorCornerIndex( c , cx , cy , cz );
			}
			else cx = cy = cz = 0;
			Stencil< Point3D< double > , 5 >& _stencil = stencils[cx][cy][cz];
			int d , off[3];
			node->depthAndOffset( d , off );
			// Set constraints from current depth
			// Gather the constraints from the vector-field at _node into the constraint stored with node
			{

				if( isInterior )
					for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
					{
						const TreeOctNode* _node = neighbors5.neighbors[x][y][z];
						if( _node )
						{
							int _idx = normalInfo.normalIndex( _node );
							if( _idx>=0 ) constraints[ node->nodeData.nodeIndex ] += Point3D< Real >::Dot( stencil.values[x][y][z] , normalInfo.normals[ _idx ] );
						}
					}
				else
					for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
					{
						const TreeOctNode* _node = neighbors5.neighbors[x][y][z];
						if( _node )
						{
							int _idx = normalInfo.normalIndex( _node );
							if( _idx>=0 )
							{
								int _d , _off[3];
								_node->depthAndOffset( _d , _off );
								constraints[ node->nodeData.nodeIndex ] += Real( GetDivergence2( integrator , d , off , _off , false , normalInfo.normals[ _idx ] ) );
								Real r = Real(GetDivergence2(integrator, d, off, _off, false, normalInfo.normals[_idx]));
							}
						}
					}
					UpdateCoarserSupportBounds( neighbors5.neighbors[2][2][2] , startX , endX , startY  , endY , startZ , endZ );
			}
			int idx = normalInfo.normalIndex( node );
			if( idx<0 ) continue;
			const Point3D< Real >& normal = normalInfo.normals[ idx ];
			if( normal[0]==0 && normal[1]==0 && normal[2]==0 ) continue;

			// Set the constraints for the parents
			if( depth>_minDepth )
			{
				neighborKey.getNeighbors( node->parent , neighbors5 );

				for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
					if( neighbors5.neighbors[x][y][z] )
					{
						TreeOctNode* _node = neighbors5.neighbors[x][y][z];
						Real c;
						if( isInterior2 )
						{
							Point3D< double >& div = _stencil.values[x][y][z];
							c = Real( div[0] * normal[0] + div[1] * normal[1] + div[2] * normal[2] );
						}
						else
						{
							int _d , _off[3];
							_node->depthAndOffset( _d , _off );
							c = Real( GetDivergence1( integrator , d , off , _off , true , normal ) );
						}
#pragma omp atomic
						_constraints[ _node->nodeData.nodeIndex ] += c;
					}
			}
		}
		MemoryUsage();
	}

	// Fine-to-coarse down-sampling of constraints
	for( int d=maxDepth-1 ; d>=(_boundaryType==0?2:0) ; d-- ) DownSample( d , _sNodes , ( ConstPointer( Real ) )_constraints + _sNodes.nodeCount[d] , _constraints+_sNodes.nodeCount[d-1] );

	// Add the accumulated constraints from all finer depths
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<_sNodes.nodeCount[maxDepth] ; i++ ) constraints[i] += _constraints[i];

	FreePointer( _constraints );


	std::vector< Point3D< Real > > coefficients( _sNodes.nodeCount[maxDepth] , zeroPoint );
	for( int d=maxDepth-1 ; d>=0 ; d-- )
	{
#pragma omp parallel for num_threads( threads )
		for( int i=_sNodes.nodeCount[d] ; i<_sNodes.nodeCount[d+1] ; i++ )
		{
			int idx = normalInfo.normalIndex( _sNodes.treeNodes[i] );
			if( idx<0 ) continue;
			coefficients[i] = normalInfo.normals[ idx ];
		}
	}

	// Coarse-to-fine up-sampling of coefficients
	for( int d=(_boundaryType==0?2:0) ; d<maxDepth ; d++ ) UpSample( d , _sNodes , ( ConstPointer( Point3D< Real > ) ) GetPointer( coefficients ) + _sNodes.nodeCount[d-1] , GetPointer( coefficients ) + _sNodes.nodeCount[d] );

	// Compute the contribution from all coarser depths
	for( int d=0 ; d<=maxDepth ; d++ )
	{
		size_t start = _sNodes.nodeCount[d] , end = _sNodes.nodeCount[d+1] , range = end - start;
		Stencil< Point3D< double > , 5 > stencils[2][2][2];
		SetDivergenceStencils( d , integrator , stencils , false );
		std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys( std::max< int >( 1 , threads ) );
		for( int i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( maxDepth );
#pragma omp parallel for num_threads( threads )
		for( int i=_sNodes.nodeCount[d] ; i<_sNodes.nodeCount[d+1] ; i++ )
		{
			typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[ omp_get_thread_num() ];
			TreeOctNode* node = _sNodes.treeNodes[i];
			int depth = node->depth();
			if( !depth ) continue;
			int startX=0 , endX=5 , startY=0 , endY=5 , startZ=0 , endZ=5;
			UpdateCoarserSupportBounds( node , startX , endX , startY  , endY , startZ , endZ );
			typename TreeOctNode::Neighbors5 neighbors5;
			neighborKey.getNeighbors( node->parent , neighbors5 );

			bool isInterior;
			{
				int d , off[3];
				node->depthAndOffset( d , off );
				int o = _boundaryType==0 ? (1<<(d-2)) : 0;
				int mn = 4+o , mx = (1<<d)-4-o;
				isInterior = ( off[0]>=mn && off[0]<mx && off[1]>=mn && off[1]<mx && off[2]>=mn && off[2]<mx );
			}
			int cx , cy , cz;
			if( d )
			{
				int c = int( node - node->parent->children );
				Cube::FactorCornerIndex( c , cx , cy , cz );
			}
			else cx = cy = cz = 0;
			Stencil< Point3D< double > , 5 >& _stencil = stencils[cx][cy][cz];

			Real constraint = Real(0);
			int d , off[3];
			node->depthAndOffset( d , off );
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
				if( neighbors5.neighbors[x][y][z] )
				{
					TreeOctNode* _node = neighbors5.neighbors[x][y][z];
					int _i = _node->nodeData.nodeIndex;
					if( isInterior )
					{
						Point3D< double >& div = _stencil.values[x][y][z];
						Point3D< Real >& normal = coefficients[_i];
						constraint += Real( div[0] * normal[0] + div[1] * normal[1] + div[2] * normal[2] );
					}
					else
					{
						int _d , _off[3];
						_node->depthAndOffset( _d , _off );
						constraint += Real( GetDivergence2( integrator , d , off , _off , true , coefficients[_i] ) );
					}
				}
				constraints[ node->nodeData.nodeIndex ] += constraint;
		}
	}
	MemoryUsage();
	return constraints;
}
template< class Real >
void Octree< Real >::refineBoundary( std::vector< int >* map ){ _sNodes.set( tree , map ); }



template< class Real >
Real Octree< Real >::getCenterValue( const typename TreeOctNode::ConstNeighborKey3& neighborKey , const TreeOctNode* node , ConstPointer( Real ) solution , ConstPointer( Real ) metSolution , const typename BSplineData< 2 >::template CenterEvaluator< 1 >& evaluator , const Stencil< double , 3 >& stencil , const Stencil< double , 3 >& pStencil , bool isInterior ) const
{
	if( node->children ) fprintf( stderr , "[WARNING] getCenterValue assumes leaf node\n" );
	Real value=0;

	int d , off[3];
	node->depthAndOffset( d , off );

	if( isInterior )
	{
		for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
		{
			const TreeOctNode* n = neighborKey.neighbors[d].neighbors[i][j][k];
			if( n ) value += solution[ n->nodeData.nodeIndex ] * Real( stencil.values[i][j][k] );
		}
		if( d>_minDepth )
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
			{
				const TreeOctNode* n = neighborKey.neighbors[d-1].neighbors[i][j][k];
				if( n ) value += metSolution[n->nodeData.nodeIndex] * Real( pStencil.values[i][j][k] );
			}
	}
	else
	{
		for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
		{
			const TreeOctNode* n = neighborKey.neighbors[d].neighbors[i][j][k];
			if( n )
			{
				int _d , _off[3];
				n->depthAndOffset( _d , _off );
				value +=
					solution[ n->nodeData.nodeIndex ] * Real(
					evaluator.value( d , off[0] , _off[0] , false , false ) * evaluator.value( d , off[1] , _off[1] , false , false ) * evaluator.value( d , off[1] , _off[1] , false , false ) );
			}
		}
		if( d>_minDepth )
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
			{
				const TreeOctNode* n = neighborKey.neighbors[d-1].neighbors[i][j][k];
				if( n )
				{
					int _d , _off[3];
					n->depthAndOffset( _d , _off );
					value +=
						solution[ n->nodeData.nodeIndex ] * Real(
						evaluator.value( d , off[0] , _off[0] , false , false ) * evaluator.value( d , off[1] , _off[1] , false , false ) * evaluator.value( d , off[1] , _off[1] , false , false ) );
				}
			}
	}
	return value;
}
template< class Real >
Real Octree< Real >::getCornerValue( const typename TreeOctNode::ConstNeighborKey3& neighborKey , const TreeOctNode* node , int corner , ConstPointer( Real ) solution , ConstPointer( Real ) metSolution , const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator , const Stencil< double , 3 >& stencil , const Stencil< double , 3 > stencils[8] , bool isInterior ) const
{
	double value = 0;
	if( _boundaryType==-1 ) value = -0.5;
	int d , off[3];
	node->depthAndOffset( d , off );

	int cx , cy , cz;
	int startX = 0 , endX = 3 , startY = 0 , endY = 3 , startZ = 0 , endZ = 3;
	Cube::FactorCornerIndex( corner , cx , cy , cz );
	{
		typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.neighbors[d];
		if( cx==0 ) endX = 2;
		else      startX = 1;
		if( cy==0 ) endY = 2;
		else      startY = 1;
		if( cz==0 ) endZ = 2;
		else      startZ = 1;
		if( isInterior )
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node=neighbors.neighbors[x][y][z];
				if( _node ) value += solution[ _node->nodeData.nodeIndex ] * stencil.values[x][y][z];
			}
		else
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if( _node )
				{
					int _d , _off[3];
					_node->depthAndOffset( _d , _off );
					value += solution[ _node->nodeData.nodeIndex ] * evaluator.value( d , off[0] , cx , _off[0] , false , false ) * evaluator.value( d , off[1] , cy , _off[1] , false , false ) * evaluator.value( d , off[2] , cz , _off[2] , false , false );
				}
			}
	}
	if( d>_minDepth )
	{
		int _corner = int( node - node->parent->children );
		int _cx , _cy , _cz;
		Cube::FactorCornerIndex( _corner , _cx , _cy , _cz );
		if( cx!=_cx ) startX = 0 , endX = 3;
		if( cy!=_cy ) startY = 0 , endY = 3;
		if( cz!=_cz ) startZ = 0 , endZ = 3;
		typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.neighbors[d-1];
		if( isInterior )
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node=neighbors.neighbors[x][y][z];
				if( _node ) value += metSolution[ _node->nodeData.nodeIndex ] * stencils[_corner].values[x][y][z];
			}
		else
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if( _node )
				{
					int _d , _off[3];
					_node->depthAndOffset( _d , _off );
					value += metSolution[ _node->nodeData.nodeIndex ] * evaluator.value( d , off[0] , cx , _off[0] , false , true ) * evaluator.value( d , off[1] , cy , _off[1] , false , true ) * evaluator.value( d , off[2] , cz , _off[2] , false , true );
				}
			}
	}
	return Real( value );
}
template< class Real >
std::pair< Real , Point3D< Real > > Octree< Real >::getCornerValueAndNormal( const typename TreeOctNode::ConstNeighborKey3& neighborKey , const TreeOctNode* node , int corner , ConstPointer( Real ) solution , ConstPointer( Real ) metSolution , const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator , const Stencil< double , 3 >& vStencil , const Stencil< double , 3 > vStencils[8] , const Stencil< Point3D< double > , 3 >& nStencil , const Stencil< Point3D< double > , 3 > nStencils[8] , bool isInterior ) const
{
	double value = 0;
	Point3D< double > normal;
	if( _boundaryType==-1 ) value = -0.5;
	int d , off[3];
	node->depthAndOffset( d , off );

	int cx , cy , cz;
	int startX = 0 , endX = 3 , startY = 0 , endY = 3 , startZ = 0 , endZ = 3;
	Cube::FactorCornerIndex( corner , cx , cy , cz );
	{
		typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.neighbors[d];
		if( cx==0 ) endX = 2;
		else      startX = 1;
		if( cy==0 ) endY = 2;
		else      startY = 1;
		if( cz==0 ) endZ = 2;
		else      startZ = 1;
		if( isInterior )
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node=neighbors.neighbors[x][y][z];
				if( _node ) value += solution[ _node->nodeData.nodeIndex ] * vStencil.values[x][y][z] , normal += nStencil.values[x][y][z] * solution[ _node->nodeData.nodeIndex ];
			}
		else
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if( _node )
				{
					int _d , _off[3];
					_node->depthAndOffset( _d , _off );
					double v [] = { evaluator.value( d , off[0] , cx , _off[0] , false , false ) , evaluator.value( d , off[1] , cy , _off[1] , false , false ) , evaluator.value( d , off[2] , cz , _off[2] , false , false ) };
					double dv[] = { evaluator.value( d , off[0] , cx , _off[0] , true  , false ) , evaluator.value( d , off[1] , cy , _off[1] , true  , false ) , evaluator.value( d , off[2] , cz , _off[2] , true  , false ) };
					value += solution[ _node->nodeData.nodeIndex ] * evaluator.value( d , off[0] , cx , _off[0] , false , false ) * evaluator.value( d , off[1] , cy , _off[1] , false , false ) * evaluator.value( d , off[2] , cz , _off[2] , false , false );
					normal += Point3D< double >( dv[0]*v[1]*v[2] , v[0]*dv[1]*v[2] , v[0]*v[1]*dv[2] ) * solution[ _node->nodeData.nodeIndex ];
				}
			}
	}
	if( d>_minDepth )
	{
		int _corner = int( node - node->parent->children );
		int _cx , _cy , _cz;
		Cube::FactorCornerIndex( _corner , _cx , _cy , _cz );
		if( cx!=_cx ) startX = 0 , endX = 3;
		if( cy!=_cy ) startY = 0 , endY = 3;
		if( cz!=_cz ) startZ = 0 , endZ = 3;
		typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.neighbors[d-1];
		if( isInterior )
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node=neighbors.neighbors[x][y][z];
				if( _node ) value += metSolution[ _node->nodeData.nodeIndex ] * vStencils[_corner].values[x][y][z] , normal += nStencils[_corner].values[x][y][z] * metSolution[ _node->nodeData.nodeIndex ];
			}
		else
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if( _node )
				{
					int _d , _off[3];
					_node->depthAndOffset( _d , _off );
					double v [] = { evaluator.value( d , off[0] , cx , _off[0] , false , true ) , evaluator.value( d , off[1] , cy , _off[1] , false , true ) , evaluator.value( d , off[2] , cz , _off[2] , false , true ) };
					double dv[] = { evaluator.value( d , off[0] , cx , _off[0] , true  , true ) , evaluator.value( d , off[1] , cy , _off[1] , true  , true ) , evaluator.value( d , off[2] , cz , _off[2] , true  , true ) };
					value += metSolution[ _node->nodeData.nodeIndex ] * evaluator.value( d , off[0] , cx , _off[0] , false , true ) * evaluator.value( d , off[1] , cy , _off[1] , false , true ) * evaluator.value( d , off[2] , cz , _off[2] , false , true );
					normal += Point3D< double >( dv[0]*v[1]*v[2] , v[0]*dv[1]*v[2] , v[0]*v[1]*dv[2] ) * metSolution[ _node->nodeData.nodeIndex ];
				}
			}
	}
	return std::pair< Real , Point3D< Real > >( Real( value ) , Point3D< Real >( normal ) );
}
template< class Real >
Point3D< Real > Octree< Real >::getCornerNormal( const typename TreeOctNode::ConstNeighbors5& neighbors5 , const typename TreeOctNode::ConstNeighbors5& pNeighbors5 , const TreeOctNode* node , int corner , ConstPointer( Real ) solution , ConstPointer( Real ) metSolution , const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator , const Stencil< Point3D< double > , 5 >& nStencil , const Stencil< Point3D< double > , 5 > nStencils[8] , bool isInterior ) const
{
	Point3D< double > normal;
	normal[0] = normal[1] = normal[2] = 0.;

	int d , off[3];
	node->depthAndOffset( d , off );

	int cx , cy , cz;
	int startX = 0 , endX = 5 , startY = 0 , endY = 5 , startZ = 0 , endZ = 5;
	Cube::FactorCornerIndex( corner , cx , cy , cz );
	{
		if( cx==0 ) endX = 4;
		else      startX = 1;
		if( cy==0 ) endY = 4;
		else      startY = 1;
		if( cz==0 ) endZ = 4;
		else      startZ = 1;
		if( isInterior )
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node=neighbors5.neighbors[x][y][z];
				if( _node ) normal += nStencil.values[x][y][z] * solution[ _node->nodeData.nodeIndex ];
			}
		else
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node = neighbors5.neighbors[x][y][z];
				if( _node )
				{
					int _d , _off[3];
					_node->depthAndOffset( _d , _off );
					double v [] = { evaluator.value( d , off[0] , cx , _off[0] , false , false ) , evaluator.value( d , off[1] , cy , _off[1] , false , false ) , evaluator.value( d , off[2] , cz , _off[2] , false , false ) };
					double dv[] = { evaluator.value( d , off[0] , cx , _off[0] , true  , false ) , evaluator.value( d , off[1] , cy , _off[1] , true  , false ) , evaluator.value( d , off[2] , cz , _off[2] , true  , false ) };
					normal += Point3D< double >( dv[0]*v[1]*v[2] , v[0]*dv[1]*v[2] , v[0]*v[1]*dv[2] ) * solution[ _node->nodeData.nodeIndex ];
				}
			}
	}
	if( d>_minDepth )
	{
		int _cx , _cy , _cz , _corner = int( node - node->parent->children );
		Cube::FactorCornerIndex( _corner , _cx , _cy , _cz );
		if( cx!=_cx ) startX = 0 , endX = 5;
		if( cy!=_cy ) startY = 0 , endY = 5;
		if( cz!=_cz ) startZ = 0 , endZ = 5;
		if( isInterior )
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node=pNeighbors5.neighbors[x][y][z];
				if( _node ) normal += nStencils[_corner].values[x][y][z] * metSolution[ _node->nodeData.nodeIndex ];
			}
		else
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node = pNeighbors5.neighbors[x][y][z];
				if( _node )
				{
					int _d , _off[3];
					_node->depthAndOffset( _d , _off );
					double v [] = { evaluator.value( d , off[0] , cx , _off[0] , false , true ) , evaluator.value( d , off[1] , cy , _off[1] , false , true ) , evaluator.value( d , off[2] , cz , _off[2] , false , true ) };
					double dv[] = { evaluator.value( d , off[0] , cx , _off[0] , true  , true ) , evaluator.value( d , off[1] , cy , _off[1] , true  , true ) , evaluator.value( d , off[2] , cz , _off[2] , true  , true ) };
					normal += Point3D< double >( dv[0]*v[1]*v[2] , v[0]*dv[1]*v[2] , v[0]*v[1]*dv[2] ) * metSolution[ _node->nodeData.nodeIndex ];
				}
			}
	}
	return Point3D< Real >( Real(normal[0]) , Real(normal[1]) , Real(normal[2]) );
}

template< class Real >
Real Octree< Real >::Evaluate( ConstPointer( Real ) coefficients , Point3D< Real > p , const BSplineData< 2 >* fData ) const
{
	Real value = Real(0);
	BSplineData< 2 > _fData;
	if( !fData ) _fData.set( tree.maxDepth() , _boundaryType ) , fData = &_fData;
	const TreeOctNode* n = tree.nextNode();
	while( n )
	{
		Point3D< Real > c;
		Real w;
		n->centerAndWidth( c , w );
		c -= p , w *= Real(1.5);
		if( fabs(c[0])>w || fabs(c[1])>w || fabs(c[2])>w )
		{
			n = tree.nextBranch( n );
			continue;
		}
		int d , off[3];
		n->depthAndOffset( d , off );
		value += (Real)
			(
			coefficients[ n->nodeData.nodeIndex ] *
			fData->baseFunctions[ BinaryNode::CenterIndex( d , off[0] ) ]( p[0] ) *
			fData->baseFunctions[ BinaryNode::CenterIndex( d , off[1] ) ]( p[1] ) *
			fData->baseFunctions[ BinaryNode::CenterIndex( d , off[2] ) ]( p[2] )
			);
		n = tree.nextNode( n );
	}
	if( _boundaryType==-1 ) value -= Real(0.5);
	return value;
}
template< class Real >
Pointer( Real ) Octree< Real >::Evaluate( ConstPointer( Real ) coefficients , int& res , Real isoValue , int depth )
{
	int maxDepth = _boundaryType==0 ? tree.maxDepth()-1 : tree.maxDepth();
	if( depth<=0 || depth>maxDepth ) depth = maxDepth;
	res = 1<<depth;
	typename BSplineData< 2 >::template ValueTables< Real > vTables = _fData.template getValueTables< Real >( _fData.VALUE_FLAG );
	Pointer( Real ) values = NewPointer< Real >( res * res * res );
	memset( values , 0 , sizeof( Real ) * res  * res * res );

	for( TreeOctNode* n=tree.nextNode() ; n ; n=tree.nextNode( n ) )
	{
		if( n->depth()>(_boundaryType==0?depth+1:depth) ) continue;
		if( n->depth()<_minDepth ) continue;
		int d , idx[3] , start[3] , end[3];
		n->depthAndOffset( d , idx );
		bool skip=false;
		for( int i=0 ; i<3 ; i++ )
		{
			// Get the index of the functions
			idx[i] = BinaryNode::CenterIndex( d , idx[i] );
			// Figure out which samples fall into the range
			vTables.setSampleSpan( idx[i] , start[i] , end[i] );
			// We only care about the odd indices
			if( !(start[i]&1) ) start[i]++;
			if( !(  end[i]&1) )   end[i]--;
			if( _boundaryType==0 )
			{
				// (start[i]-1)>>1 >=   res/2 
				// (  end[i]-1)<<1 <  3*res/2
				start[i] = std::max< int >( start[i] ,   res+1 );
				end  [i] = std::min< int >( end  [i] , 3*res-1 );
			}
		}
		if( skip ) continue;
		Real coefficient = coefficients[ n->nodeData.nodeIndex ];
		for( int x=start[0] ; x<=end[0] ; x+=2 )
			for( int y=start[1] ; y<=end[1] ; y+=2 )
				for( int z=start[2] ; z<=end[2] ; z+=2 )
				{
					int xx = (x-1)>>1 , yy=(y-1)>>1 , zz = (z-1)>>1;
					if( _boundaryType==0 ) xx -= res/2 , yy -= res/2 , zz -= res/2;
					values[ zz*res*res + yy*res + xx ] +=
						coefficient *
						vTables.valueTable[ idx[0] + x*vTables.functionCount ] *
						vTables.valueTable[ idx[1] + y*vTables.functionCount ] *
						vTables.valueTable[ idx[2] + z*vTables.functionCount ];
				}
	}
	if( _boundaryType==-1 ) for( int i=0 ; i<res*res*res ; i++ ) values[i] -= Real(0.5);
	for( int i=0 ; i<res*res*res ; i++ ) values[i] -= isoValue;

	return values;
}

//calculating the gradient of the implicit function.
template< class Real >
void Octree< Real >::ImplicitGrad(ConstPointer(Real) coefficients, Point3D< Real > p, Real* grad_p, BSplineData< 2 >* fData) const
{
	Real value1 = Real(0);
	Real value2 = Real(0);
	Real value3 = Real(0);
	BSplineData< 2 > _fData;
	if (!fData) _fData.set(tree.maxDepth(), _boundaryType), fData = &_fData;
	const TreeOctNode* n = tree.nextNode();
	while (n)
	{
		Point3D< Real > c;
		Real w;
		n->centerAndWidth(c, w);
		c -= p, w *= Real(1.5);
		if (fabs(c[0]) > w || fabs(c[1]) > w || fabs(c[2]) > w)
		{
			n = tree.nextBranch(n);
			continue;
		}
		int d, off[3];
		n->depthAndOffset(d, off);
		value1 += (Real)
			(
				coefficients[n->nodeData.nodeIndex] *
				(fData->baseFunctions[BinaryNode::CenterIndex(d, off[0])].derivative())(p[0]) *
				fData->baseFunctions[BinaryNode::CenterIndex(d, off[1])](p[1]) *
				fData->baseFunctions[BinaryNode::CenterIndex(d, off[2])](p[2])
				);
		value2 += (Real)
			(
				coefficients[n->nodeData.nodeIndex] *
				fData->baseFunctions[BinaryNode::CenterIndex(d, off[0])](p[0]) *
				(fData->baseFunctions[BinaryNode::CenterIndex(d, off[1])].derivative())(p[1]) *
				fData->baseFunctions[BinaryNode::CenterIndex(d, off[2])](p[2])
				);
		value3 += (Real)
			(
				coefficients[n->nodeData.nodeIndex] *
				fData->baseFunctions[BinaryNode::CenterIndex(d, off[0])](p[0]) *
				fData->baseFunctions[BinaryNode::CenterIndex(d, off[1])](p[1]) *
				(fData->baseFunctions[BinaryNode::CenterIndex(d, off[2])].derivative())(p[2])
				);
		n = tree.nextNode(n);
	}
	grad_p[0] = value1;
	grad_p[1] = value2;
	grad_p[2] = value3;
	return;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Octree2. Containing the main implementation of the coarse point cloud orientation module.
////////////
template< class Real > double Octree2< Real >::maxMemoryUsage = 0;

template< class Real >
double Octree2< Real >::MemoryUsage(void)
{
	double mem = double(MemoryInfo::Usage()) / (1 << 20);
	if (mem > maxMemoryUsage) maxMemoryUsage = mem;
	return mem;
}

template< class Real >
Octree2< Real >::Octree2(void)
{
	num_point = 0;
	threads = 1;
	_normalSmooth = 0;
	_constrainValues = false;
}

template< class Real >
bool Octree2< Real >::_IsInset(const TreeOctNode* node)
{
	int d, off[3];
	node->depthAndOffset(d, off);
	int res = 1 << d, o = 1 << (d - 2);
	return (off[0] >= o && off[0] < res - o && off[1] >= o && off[1] < res - o && off[2] >= o && off[2] < res - o);
}
template< class Real >
bool Octree2< Real >::_IsInsetSupported(const TreeOctNode* node)
{
	int d, off[3];
	node->depthAndOffset(d, off);
	int res = 1 << d, o = (1 << (d - 2)) - 1;
	return (off[0] >= o && off[0] < res - o && off[1] >= o && off[1] < res - o && off[2] >= o && off[2] < res - o);
}
template< class Real >
int Octree2< Real >::SplatOrientedPoint(ConstPointer(Real) kernelDensityWeights, TreeOctNode* node, const Point3D<Real>& position, Real& decoeff, SmoothCoeff& smoothcoeff, int cnt, typename TreeOctNode::NeighborKey3& neighborKey)
{
	double x, dxdy, dxdydz, dx[DIMENSION][SPLAT_ORDER + 1];
	double width;
	int off[3];
	typename TreeOctNode::Neighbors3& neighbors = neighborKey.setNeighbors(node);
	Point3D<Real> center;
	Real w;
	node->centerAndWidth(center, w);
	width = w;
	for (int i = 0; i < 3; i++)
	{
#if SPLAT_ORDER==2
		off[i] = 0;
		x = (center[i] - position[i] - width) / width;
		dx[i][0] = 1.125 + 1.500 * x + 0.500 * x * x;
		x = (center[i] - position[i]) / width;
		dx[i][1] = 0.750 - x * x;

		dx[i][2] = 1. - dx[i][1] - dx[i][0];
#elif SPLAT_ORDER==1
		x = (position[i] - center[i]) / width;
		if (x < 0)
		{
			off[i] = 0;
			dx[i][0] = -x;
		}
		else
		{
			off[i] = 1;
			dx[i][0] = 1. - x;
		}
		dx[i][1] = 1. - dx[i][0];
#elif SPLAT_ORDER==0
		off[i] = 1;
		dx[i][0] = 1.;
#else
#     error Splat order not supported
#endif // SPLAT_ORDER
	}
	for (int i = off[0]; i <= off[0] + SPLAT_ORDER; i++) for (int j = off[1]; j <= off[1] + SPLAT_ORDER; j++)
	{
		dxdy = dx[0][i] * dx[1][j];
		for (int k = off[2]; k <= off[2] + SPLAT_ORDER; k++)
			if (neighbors.neighbors[i][j][k])
			{
				dxdydz = dxdy * dx[2][k];
				TreeOctNode* _node = neighbors.neighbors[i][j][k];
				if (smoothcoeff.scoeffIndices.size() < TreeNodeData::NodeCount) smoothcoeff.scoeffIndices.resize(TreeNodeData::NodeCount, -1);
				int idx = smoothcoeff.scoeffIndex(_node);
				if (idx < 0)
				{
					idx = smoothcoeff.scoeffIndices[_node->nodeData.nodeIndex] = (int)smoothcoeff.smoothcoeff.size();
					smoothcoeff.smoothcoeff.push_back(std::vector< MatrixEntry <Real> > ());
					smoothcoeff.weight.push_back(Real(0));
					smoothcoeff.boundary.push_back(std::vector<bool>(3, false));
					/*
					Point3D<Real> n;
					n[0] = n[1] = n[2] = 0;
					idx = normalInfo.normalIndices[_node->nodeData.nodeIndex] = (int)normalInfo.normals.size();
					normalInfo.normals.push_back(n);
					*/
				}
				//smoothcoeff.smoothcoeff[idx].push_back(MatrixEntry< Real >(cnt, Real(dxdydz)));
				//smoothcoeff.weight[idx] += Real(dxdydz);
				Real thre = 1e-6;

				if (Real(abs(dxdydz * decoeff)) > thre && !isnan(decoeff))
				{
					std::vector<MatrixEntry <Real>> &l = smoothcoeff.smoothcoeff[idx];
					if (l.size() == 0)
					{
						l.push_back(MatrixEntry< Real >(cnt, Real(dxdydz * decoeff)));
					}
					else
					{
						if (l[l.size() - 1].N == cnt)
						{
							l[l.size() - 1].Value += Real(dxdydz * decoeff);
						}
						else
						{
							l.push_back(MatrixEntry< Real >(cnt, Real(dxdydz * decoeff)));
						}
					}
					smoothcoeff.weight[idx] += Real(dxdydz * decoeff);
				}
			}
	}
	return 0;
}
template< class Real >
Real Octree2< Real >::SplatOrientedPoint(ConstPointer(Real) kernelDensityWeights, const Point3D<Real>& position, SmoothCoeff& smoothcoeff, int cnt, typename TreeOctNode::NeighborKey3& neighborKey, int splatDepth, Real samplesPerNode, int minDepth, int maxDepth)
{
	double dx;
	Real dencoeff;
	Point3D<Real> n;
	TreeOctNode* temp;
	double width;
	Point3D< Real > myCenter;
	Real myWidth;
	myCenter[0] = myCenter[1] = myCenter[2] = Real(0.5);
	myWidth = Real(1.0);

	temp = &tree;
	while (temp->depth() < splatDepth)
	{
		if (!temp->children)
		{
			fprintf(stderr, "Octree::SplatOrientedPoint error\n");
			return -1;
		}
		int cIndex = TreeOctNode::CornerIndex(myCenter, position);
		temp = &temp->children[cIndex];
		myWidth /= 2;
		if (cIndex & 1) myCenter[0] += myWidth / 2;
		else		 myCenter[0] -= myWidth / 2;
		if (cIndex & 2) myCenter[1] += myWidth / 2;
		else		 myCenter[1] -= myWidth / 2;
		if (cIndex & 4) myCenter[2] += myWidth / 2;
		else		 myCenter[2] -= myWidth / 2;
	}
	Real weight, depth;
	GetSampleDepthAndWeight(kernelDensityWeights, temp, position, neighborKey, samplesPerNode, depth, weight);

	if (depth < minDepth) depth = Real(minDepth);
	if (depth > maxDepth) depth = Real(maxDepth);
	int topDepth = int(ceil(depth));

	dx = 1.0 - (topDepth - depth);
	if (topDepth <= minDepth)
	{
		topDepth = minDepth;
		dx = 1;
	}
	else if (topDepth > maxDepth)
	{
		topDepth = maxDepth;
		dx = 1;
	}
	while (temp->depth() > topDepth) temp = temp->parent;
	while (temp->depth() < topDepth)
	{
		if (!temp->children) temp->initChildren();
		int cIndex = TreeOctNode::CornerIndex(myCenter, position);
		temp = &temp->children[cIndex];
		myWidth /= 2;
		if (cIndex & 1) myCenter[0] += myWidth / 2;
		else		 myCenter[0] -= myWidth / 2;
		if (cIndex & 2) myCenter[1] += myWidth / 2;
		else		 myCenter[1] -= myWidth / 2;
		if (cIndex & 4) myCenter[2] += myWidth / 2;
		else		 myCenter[2] -= myWidth / 2;
	}
	width = 1.0 / (1 << temp->depth());
	dencoeff = weight / Real(pow(width, 3)) * Real(dx);
	//n = normal * weight / Real(pow(width, 3)) * Real(dx);
	SplatOrientedPoint(kernelDensityWeights, temp, position, dencoeff, smoothcoeff, cnt, neighborKey);
	if (fabs(1.0 - dx) > EPSILON)
	{
		dx = Real(1.0 - dx);
		temp = temp->parent;
		width = 1.0 / (1 << temp->depth());

		//n = normal * weight / Real(pow(width, 3)) * Real(dx);
		dencoeff = weight / Real(pow(width, 3)) * Real(dx);
		SplatOrientedPoint(kernelDensityWeights, temp, position, dencoeff, smoothcoeff, cnt, neighborKey);
	}
	return weight;
}

template< class Real >
void Octree2< Real >::GetSampleDepthAndWeight(ConstPointer(Real) kernelDensityWeights, const TreeOctNode* node, const Point3D<Real>& position, typename TreeOctNode::ConstNeighborKey3& neighborKey, Real samplesPerNode, Real& depth, Real& weight)
{
	const TreeOctNode* temp = node;
	weight = Real(1.0) / GetSampleWeight(kernelDensityWeights, temp, position, neighborKey);
	if (weight >= samplesPerNode) depth = Real(temp->depth() + log(weight / samplesPerNode) / log(double(1 << (DIMENSION - 1))));
	else
	{
		Real oldWeight, newWeight;
		oldWeight = newWeight = weight;
		while (newWeight < samplesPerNode && temp->parent)
		{
			temp = temp->parent;
			oldWeight = newWeight;
			newWeight = Real(1.0) / GetSampleWeight(kernelDensityWeights, temp, position, neighborKey);
		}
		depth = Real(temp->depth() + log(newWeight / samplesPerNode) / log(newWeight / oldWeight));
	}
	weight = Real(pow(double(1 << (DIMENSION - 1)), -double(depth)));
}
template< class Real >
void Octree2< Real >::GetSampleDepthAndWeight(ConstPointer(Real) kernelDensityWeights, TreeOctNode* node, const Point3D<Real>& position, typename TreeOctNode::NeighborKey3& neighborKey, Real samplesPerNode, Real& depth, Real& weight)
{
	TreeOctNode* temp = node;
	weight = Real(1.0) / GetSampleWeight(kernelDensityWeights, temp, position, neighborKey);
	if (weight >= samplesPerNode) depth = Real(temp->depth() + log(weight / samplesPerNode) / log(double(1 << (DIMENSION - 1))));
	else
	{
		Real oldWeight, newWeight;
		oldWeight = newWeight = weight;
		while (newWeight < samplesPerNode && temp->parent)
		{
			temp = temp->parent;
			oldWeight = newWeight;
			newWeight = Real(1.0) / GetSampleWeight(kernelDensityWeights, temp, position, neighborKey);
		}
		depth = Real(temp->depth() + log(newWeight / samplesPerNode) / log(newWeight / oldWeight));
	}
	weight = Real(pow(double(1 << (DIMENSION - 1)), -double(depth)));
}
template< class Real >
Real Octree2< Real >::GetSampleWeight(ConstPointer(Real) kernelDensityWeights, const Point3D<Real>& position, typename TreeOctNode::NeighborKey3& neighborKey, int splatDepth)
{
	Point3D< Real > myCenter;
	Real myWidth;
	myCenter[0] = myCenter[1] = myCenter[2] = Real(0.5);
	myWidth = Real(1.0);

	TreeOctNode* temp = &tree;
	int d = 0;
	while (d < splatDepth)
	{
		if (!temp->children)
		{
			fprintf(stderr, "Octree::SplatOrientedPoint error\n");
			return -1;
		}
		int cIndex = TreeOctNode::CornerIndex(myCenter, position);
		temp = &temp->children[cIndex];
		myWidth /= 2;
		if (cIndex & 1) myCenter[0] += myWidth / 2;
		else 		   myCenter[0] -= myWidth / 2;
		if (cIndex & 2) myCenter[1] += myWidth / 2;
		else 		   myCenter[1] -= myWidth / 2;
		if (cIndex & 4) myCenter[2] += myWidth / 2;
		else 		   myCenter[2] -= myWidth / 2;
		d++;
	}
	return GetSampleWeight(kernelDensityWeights, temp, position, neighborKey);
}
//GetSampleWeight是插值附近node的kernelDensityWeights，并返回其倒数
template< class Real >
Real Octree2< Real >::GetSampleWeight(ConstPointer(Real) kernelDensityWeights, TreeOctNode* node, const Point3D<Real>& position, typename TreeOctNode::NeighborKey3& neighborKey)
{
	Real weight = 0;
	double x, dxdy, dx[DIMENSION][3];
	double width;
	typename TreeOctNode::Neighbors3& neighbors = neighborKey.setNeighbors(node);
	Point3D<Real> center;
	Real w;
	node->centerAndWidth(center, w);
	width = w;

	for (int i = 0; i < DIMENSION; i++)
	{
		x = (center[i] - position[i] - width) / width;
		dx[i][0] = 1.125 + 1.500 * x + 0.500 * x * x;
		x = (center[i] - position[i]) / width;
		dx[i][1] = 0.750 - x * x;

		dx[i][2] = 1.0 - dx[i][1] - dx[i][0];
	}

	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
	{
		dxdy = dx[0][i] * dx[1][j];
		for (int k = 0; k < 3; k++) if (neighbors.neighbors[i][j][k])
			weight += Real(dxdy * dx[2][k] * kernelDensityWeights[neighbors.neighbors[i][j][k]->nodeData.nodeIndex]);
	}
	return Real(1.0 / weight);
}
template< class Real >
Real Octree2< Real >::GetSampleWeight(ConstPointer(Real) kernelDensityWeights, const TreeOctNode* node, const Point3D<Real>& position, typename TreeOctNode::ConstNeighborKey3& neighborKey)
{
	Real weight = 0;
	double x, dxdy, dx[DIMENSION][3];
	double width;
	typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.getNeighbors(node);
	Point3D<Real> center;
	Real w;
	node->centerAndWidth(center, w);
	width = w;
	//i->dimension j->contribution of two neighbours and itself
	for (int i = 0; i < DIMENSION; i++)
	{
		x = (center[i] - position[i] - width) / width;
		dx[i][0] = 1.125 + 1.500 * x + 0.500 * x * x;
		x = (center[i] - position[i]) / width;
		dx[i][1] = 0.750 - x * x;

		dx[i][2] = 1.0 - dx[i][1] - dx[i][0];
	}

	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
	{
		dxdy = dx[0][i] * dx[1][j];
		for (int k = 0; k < 3; k++) if (neighbors.neighbors[i][j][k])
			weight += Real(dxdy * dx[2][k] * kernelDensityWeights[neighbors.neighbors[i][j][k]->nodeData.nodeIndex]);
	}
	return Real(1.0 / weight);
}
template< class Real >
int Octree2< Real >::UpdateWeightContribution(std::vector< Real >& kernelDensityWeights, TreeOctNode* node, const Point3D<Real>& position, typename TreeOctNode::NeighborKey3& neighborKey, Real weight)
{
	typename TreeOctNode::Neighbors3& neighbors = neighborKey.setNeighbors(node);
	if (kernelDensityWeights.size() < TreeNodeData::NodeCount) kernelDensityWeights.resize(TreeNodeData::NodeCount, 0);
	double x, dxdy, dx[DIMENSION][3], width;
	Point3D< Real > center;
	Real w;
	node->centerAndWidth(center, w);
	width = w;
	const double SAMPLE_SCALE = 1. / (0.125 * 0.125 + 0.75 * 0.75 + 0.125 * 0.125);

	for (int i = 0; i < DIMENSION; i++)
	{
		x = (center[i] - position[i] - width) / width;
		dx[i][0] = 1.125 + 1.500 * x + 0.500 * x * x;
		dx[i][1] = -0.25 - 2. * x - x * x;
		dx[i][2] = 1. - dx[i][1] - dx[i][0];
		// Note that we are splatting along a co-dimension one manifold, so uniform point samples
		// do not generate a unit sample weight.
		dx[i][0] *= SAMPLE_SCALE;
	}
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
	{
		dxdy = dx[0][i] * dx[1][j] * weight;
		TreeOctNode** _neighbors = neighbors.neighbors[i][j];
		for (int k = 0; k < 3; k++) if (_neighbors[k]) kernelDensityWeights[_neighbors[k]->nodeData.nodeIndex] += Real(dxdy * dx[2][k]);
	}
	return 0;
}
template< class Real >
bool Octree2< Real >::_InBounds(Point3D< Real > p) const
{
	if (_boundaryType == 0) { if (p[0]<Real(0.25) || p[0]>Real(0.75) || p[1]<Real(0.25) || p[1]>Real(0.75) || p[2]<Real(0.25) || p[2]>Real(0.75)) return false; }
	else { if (p[0]<Real(0.00) || p[0]>Real(1.00) || p[1]<Real(0.00) || p[1]>Real(1.00) || p[2]<Real(0.00) || p[2]>Real(1.00)) return false; }
	return true;
}
template< class Real >
template< class PointReal >
int Octree2< Real >::SetTree(PointStream< PointReal >* pointStream, int minDepth, int maxDepth, int fullDepth,
	int splatDepth, Real samplesPerNode, Real scaleFactor, int adaptiveExponent,
	PointInfo& pointInfo, NormalInfo& normalInfo, SmoothCoeff &smoothcoeff, std::vector< Real >& kernelDensityWeights, std::vector< Real >& centerWeights,
	int boundaryType, XForm4x4< Real > xForm, bool makeComplete)
{
	if (splatDepth < 0) splatDepth = 0;

	_boundaryType = boundaryType;
	if (_boundaryType < 0) _boundaryType = -1;
	else if (_boundaryType > 0) _boundaryType = 1;
	_samplesPerNode = samplesPerNode;
	_splatDepth = splatDepth;
	XForm3x3< Real > xFormN;
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) xFormN(i, j) = xForm(i, j);
	xFormN = xFormN.transpose().inverse();
	minDepth = std::min< int >(minDepth, maxDepth);	// minDepth <= maxDepth
	fullDepth = std::max< int >(minDepth, std::min< int >(fullDepth, maxDepth));	// minDepth <= fullDepth <= maxDepth
	// If _boundaryType==0, points are scaled to be in the [0.25,0.75]^3 cube so all depths have to be offset by
	// and the minDepth has to be 2.
	if (_boundaryType == 0)
	{
		minDepth++, maxDepth++, fullDepth++;
		if (splatDepth) splatDepth++;
		minDepth = std::max< int >(minDepth, 2);
	}
	// Otherwise the points are in the [0,1]^3 cube.
	// However, for Neumann constraints, the function at depth 0 is constant so the system matrix is zero if there
	// is no screening.
	_fData.set(maxDepth, _boundaryType);

	_minDepth = minDepth;
	_fullDepth = fullDepth;
	double pointWeightSum = 0;
	Point3D< Real > min, max, myCenter;
	Real myWidth;
	int i, cnt = 0;
	TreeOctNode* temp;

	typename TreeOctNode::NeighborKey3 neighborKey;
	neighborKey.set(maxDepth);
	tree.setFullDepth(_fullDepth);
	
	// Read through once to get the center
	{
		double t = Time();
		Point3D< Real > p, n;
		Point3D< PointReal > _p, _n;
		while (pointStream->nextPoint(_p, _n))
		{
			p = xForm * Point3D< Real >(_p);
			for (i = 0; i < DIMENSION; i++)
			{
				if (!cnt || p[i] < min[i]) min[i] = p[i];
				if (!cnt || p[i] > max[i]) max[i] = p[i];
			}
			cnt++;
		}

		if (_boundaryType == 0) _scale = std::max< Real >(max[0] - min[0], std::max< Real >(max[1] - min[1], max[2] - min[2])) * 2;
		else         _scale = std::max< Real >(max[0] - min[0], std::max< Real >(max[1] - min[1], max[2] - min[2]));
		_center = (max + min) / 2;
	}
	_scaleFactor = scaleFactor;
	_scale *= scaleFactor;
	for (i = 0; i < DIMENSION; i++)
	{
		_center[i] -= _scale / 2;
		//printf("%f\n", _center[i]);
	}
	if (splatDepth > 0)
	{
		double t = Time();
		cnt = 0;
		pointStream->reset();
		Point3D< Real > p, n;
		Point3D< PointReal > _p, _n;
		while (pointStream->nextPoint(_p, _n))
		{
			p = xForm * Point3D< Real >(_p), n = xFormN * Point3D< Real >(_n);;
			p = (p - _center) / _scale;
			if (!_InBounds(p)) continue;
			myCenter = Point3D< Real >(Real(0.5), Real(0.5), Real(0.5));
			myWidth = Real(1.0);
			Real weight = Real(1.);
			temp = &tree;
			int d = 0;

			while (d < splatDepth)
			{
				UpdateWeightContribution(kernelDensityWeights, temp, p, neighborKey, weight);
				if (!temp->children) temp->initChildren();
				int cIndex = TreeOctNode::CornerIndex(myCenter, p);
				temp = temp->children + cIndex;
				myWidth /= 2;
				if (cIndex & 1) myCenter[0] += myWidth / 2;
				else           myCenter[0] -= myWidth / 2;
				if (cIndex & 2) myCenter[1] += myWidth / 2;
				else           myCenter[1] -= myWidth / 2;
				if (cIndex & 4) myCenter[2] += myWidth / 2;
				else           myCenter[2] -= myWidth / 2;
				d++;
			}
			UpdateWeightContribution(kernelDensityWeights, temp, p, neighborKey, weight);
			cnt++;
		}
	}
	kernelDensityWeights.resize(TreeNodeData::NodeCount, 0);

	std::vector< _PointData >& points = pointInfo.points;
	num_point = cnt;
	cnt = 0;
	pointStream->reset();
	Point3D< Real > p, n;
	Point3D< PointReal > _p, _n;
	while (pointStream->nextPoint(_p, _n))
	{
		p = xForm * Point3D< Real >(_p), n = xFormN * Point3D< Real >(_n);
		p = (p - _center) / _scale;
		n *= Real(-1.);
		all_points.push_back(p);
		all_normals.push_back(n);
		//printf("%d\n", all_points.size());
		if (!_InBounds(p)) continue;
		myCenter = Point3D< Real >(Real(0.5), Real(0.5), Real(0.5));
		myWidth = Real(1.0);


		Real pointWeight = Real(1.f);
		//pointWeight equals area P 
		if (samplesPerNode > 0 && splatDepth) pointWeight = SplatOrientedPoint(GetPointer(kernelDensityWeights), p, smoothcoeff, cnt, neighborKey, splatDepth, samplesPerNode, _minDepth, maxDepth);
		else
		{
			temp = &tree;
			int d = 0;
			if (splatDepth)
			{
				while (d < splatDepth)
				{
					int cIndex = TreeOctNode::CornerIndex(myCenter, p);
					temp = &temp->children[cIndex];
					myWidth /= 2;
					if (cIndex & 1) myCenter[0] += myWidth / 2;
					else		 myCenter[0] -= myWidth / 2;
					if (cIndex & 2) myCenter[1] += myWidth / 2;
					else		 myCenter[1] -= myWidth / 2;
					if (cIndex & 4) myCenter[2] += myWidth / 2;
					else		 myCenter[2] -= myWidth / 2;
					d++;
				}
				pointWeight = GetSampleWeight(GetPointer(kernelDensityWeights), temp, p, neighborKey);
			}
			//for (i = 0; i < DIMENSION; i++) n[i] *= pointWeight;
			while (d < maxDepth)
			{
				if (!temp->children) temp->initChildren();
				int cIndex = TreeOctNode::CornerIndex(myCenter, p);
				temp = &temp->children[cIndex];
				myWidth /= 2;
				if (cIndex & 1) myCenter[0] += myWidth / 2;
				else		 myCenter[0] -= myWidth / 2;
				if (cIndex & 2) myCenter[1] += myWidth / 2;
				else		 myCenter[1] -= myWidth / 2;
				if (cIndex & 4) myCenter[2] += myWidth / 2;
				else		 myCenter[2] -= myWidth / 2;
				d++;
			}
			SplatOrientedPoint(GetPointer(kernelDensityWeights), temp, p, pointWeight, smoothcoeff, cnt, neighborKey);
		}
		pointWeightSum += pointWeight;
		
		if (1)
		{
			Real pointScreeningWeight = Real(1.f);
			int d = 0;
			TreeOctNode* temp = &tree;
			myCenter = Point3D< Real >(Real(0.5), Real(0.5), Real(0.5));
			myWidth = Real(1.0);
			while (1)
			{
				if (pointInfo.pointIndices.size() < TreeNodeData::NodeCount) pointInfo.pointIndices.resize(TreeNodeData::NodeCount, -1);
				int idx = pointInfo.pointIndex(temp);

				if (idx == -1)
				{
					idx = (int)points.size();
					points.push_back(_PointData(p * pointScreeningWeight, pointScreeningWeight));
					pointInfo.pointIndices[temp->nodeData.nodeIndex] = idx;
				}
				else
				{
					points[idx].weight += pointScreeningWeight;
					points[idx].position += p * pointScreeningWeight;
				}

				int cIndex = TreeOctNode::CornerIndex(myCenter, p);
				if (!temp->children)
				{
					temp->nodeData.leaf_to_pidx.push_back(cnt);
					break;
				}
				temp = &temp->children[cIndex];
				myWidth /= 2;
				if (cIndex & 1) myCenter[0] += myWidth / 2;
				else		   myCenter[0] -= myWidth / 2;
				if (cIndex & 2) myCenter[1] += myWidth / 2;
				else		   myCenter[1] -= myWidth / 2;
				if (cIndex & 4) myCenter[2] += myWidth / 2;
				else		   myCenter[2] -= myWidth / 2;
				d++;
			}
		}
		cnt++;
	}
	//printf("point_size:%d\n", (int)points.size());
	Real constraintWeight = 4.0;
	if (_boundaryType == 0) pointWeightSum *= Real(4.);
	constraintWeight *= Real(pointWeightSum);
	constraintWeight /= cnt;

	MemoryUsage();
	// Set the average position and scale the weights
	for (TreeOctNode* node = tree.nextNode(); node; node = tree.nextNode(node))
		if (pointInfo.pointIndex(node) != -1)
		{
			int idx = pointInfo.pointIndex(node);
			points[idx].position /= points[idx].weight;
			int e = (_boundaryType == 0 ? node->depth() - 1 : node->depth()) * adaptiveExponent - (_boundaryType == 0 ? maxDepth - 1 : maxDepth) * (adaptiveExponent - 1);
			if (e < 0) points[idx].weight /= Real(1 << (-e));
			else      points[idx].weight *= Real(1 << e);
			points[idx].weight *= Real(constraintWeight);
		}

#if FORCE_NEUMANN_FIELD
	if (_boundaryType == 1)
		for (TreeOctNode* node = tree.nextNode(); node; node = tree.nextNode(node))
		{
			int d, off[3], res;
			node->depthAndOffset(d, off);
			res = 1 << d;
			int idx = smoothcoeff.scoeffIndex(node);
			if (idx < 0) continue;
			//std::vector< MatrixEntry <Real> >& coeff = smoothcoeff.smoothcoeff[idx];
			for (int dim = 0; dim < 3; dim++)
			{
				if (off[dim] == 0 || off[dim] == res - 1)
				{
					smoothcoeff.boundary[idx][dim] = true;
				}
			}
		}
#endif // FORCE_NEUMANN_FIELD

	centerWeights.resize(tree.nodes(), 0);
	kernelDensityWeights.resize(tree.nodes(), 0);
	// Set the point weights for evaluating the iso-value
	for (TreeOctNode* node = tree.nextNode(); node; node = tree.nextNode(node))
	{
		int idx = smoothcoeff.scoeffIndex(node);
		if (idx < 0) centerWeights[node->nodeData.nodeIndex] = 0;
		else        centerWeights[node->nodeData.nodeIndex] = smoothcoeff.weight[idx];
	}
	MemoryUsage();
	{
		std::vector< int > indexMap;
		if (makeComplete) MakeComplete(&indexMap);
		else ClipTree(smoothcoeff), Finalize(&indexMap);
		
		{
			std::vector< int > temp = pointInfo.pointIndices;
			pointInfo.pointIndices.resize(indexMap.size());
			for (int i = 0; i < indexMap.size(); i++)
				if (indexMap[i] < temp.size()) pointInfo.pointIndices[i] = temp[indexMap[i]];
				else                          pointInfo.pointIndices[i] = -1;
		}
		
		{
			std::vector< int > temp = smoothcoeff.scoeffIndices;
			smoothcoeff.scoeffIndices.resize(indexMap.size());
			for (int i = 0; i < indexMap.size(); i++)
				if (indexMap[i] < temp.size()) smoothcoeff.scoeffIndices[i] = temp[indexMap[i]];
				else                          smoothcoeff.scoeffIndices[i] = -1;
		}
		{
			std::vector< Real > temp = centerWeights;
			centerWeights.resize(indexMap.size());
			for (int i = 0; i < indexMap.size(); i++)
				if (indexMap[i] < temp.size()) centerWeights[i] = temp[indexMap[i]];
				else                          centerWeights[i] = (Real)0;
		}
		{
			std::vector< Real > temp = kernelDensityWeights;
			kernelDensityWeights.resize(indexMap.size());
			for (int i = 0; i < indexMap.size(); i++)
				if (indexMap[i] < temp.size()) kernelDensityWeights[i] = temp[indexMap[i]];
				else                          kernelDensityWeights[i] = (Real)0;
		}
	}
	return cnt;
}
template< class Real >
void Octree2< Real >::MakeComplete(std::vector< int >* map)
{
	tree.setFullDepth(tree.maxDepth());
	refineBoundary(map);
	MemoryUsage();
}
template< class Real >
void Octree2< Real >::ClipTree(const SmoothCoeff& smoothcoeff)
{
	int maxDepth = tree.maxDepth();
	for (TreeOctNode* temp = tree.nextNode(); temp; temp = tree.nextNode(temp))
		if (temp->children && temp->depth() >= _fullDepth)
		{
			int hasNormals = 0;
			for (int i = 0; i < Cube::CORNERS && !hasNormals; i++) hasNormals = HasNormals(&temp->children[i], smoothcoeff);
			if (!hasNormals) temp->children = NULL;
		}
	MemoryUsage();
}

template< class Real >
void Octree2< Real >::Finalize(std::vector< int >* map)
{
	int maxDepth = tree.maxDepth();
	typename TreeOctNode::NeighborKey3 neighborKey;
	neighborKey.set(maxDepth);
	for (int d = maxDepth; d > 1; d--)
		for (TreeOctNode* node = tree.nextNode(); node; node = tree.nextNode(node)) if (node->depth() == d)
		{
			typename TreeOctNode::Neighbors3& neighbors = neighborKey.setNeighbors(node->parent->parent);
			for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++)
				if (neighbors.neighbors[i][j][k] && !neighbors.neighbors[i][j][k]->children)
					neighbors.neighbors[i][j][k]->initChildren();
		}
	refineBoundary(map);
}
template< class Real >
double Octree2< Real >::GetLaplacian(const typename BSplineData< 2 >::Integrator& integrator, int d, const int off1[], const int off2[], bool childParent) const
{
	double vv[] =
	{
		integrator.dot(d , off1[0] , off2[0] , false , false , childParent) ,
		integrator.dot(d , off1[1] , off2[1] , false , false , childParent) ,
		integrator.dot(d , off1[2] , off2[2] , false , false , childParent)
	};
	double dd[] =
	{
		integrator.dot(d , off1[0] , off2[0] , true , true , childParent) ,
		integrator.dot(d , off1[1] , off2[1] , true , true , childParent) ,
		integrator.dot(d , off1[2] , off2[2] , true , true , childParent)
	};
	return dd[0] * vv[1] * vv[2] + vv[0] * dd[1] * vv[2] + vv[0] * vv[1] * dd[2];
}
template< class Real >
double Octree2< Real >::GetDivergence1(const typename BSplineData< 2 >::Integrator& integrator, int d, const int off1[], const int off2[], bool childParent, const Point3D< Real >& normal1) const
{
	return Point3D< double >::Dot(GetDivergence1(integrator, d, off1, off2, childParent), normal1);
}
template< class Real >
double Octree2< Real >::GetDivergence2(const typename BSplineData< 2 >::Integrator& integrator, int d, const int off1[], const int off2[], bool childParent, const Point3D< Real >& normal2) const
{
	return Point3D< double >::Dot(GetDivergence2(integrator, d, off1, off2, childParent), normal2);
}
template< class Real >
Point3D< double > Octree2< Real >::GetDivergence1(const typename BSplineData< 2 >::Integrator& integrator, int d, const int off1[], const int off2[], bool childParent) const
{
	//B(p) = B(x)B(y)B(z), grad(B(p)) = grad(B(x))B(y)B(z) + ... + ...
	double vv[] =
	{
		integrator.dot(d , off1[0] , off2[0] , false , false , childParent) ,
		integrator.dot(d , off1[1] , off2[1] , false , false , childParent) ,
		integrator.dot(d , off1[2] , off2[2] , false , false , childParent)
	};
#if GRADIENT_DOMAIN_SOLUTION
	// Take the dot-product of the vector-field with the gradient of the basis function
	double vd[] =
	{
		integrator.dot(d , off1[0] , off2[0] , false , true , childParent) ,
		integrator.dot(d , off1[1] , off2[1] , false , true , childParent) ,
		integrator.dot(d , off1[2] , off2[2] , false , true , childParent)
	};
	return  Point3D< double >(vd[0] * vv[1] * vv[2], vv[0] * vd[1] * vv[2], vv[0] * vv[1] * vd[2]);
#else // !GRADIENT_DOMAIN_SOLUTION
	// Take the dot-product of the divergence of the vector-field with the basis function
	double dv[] =
	{
		integrator.dot(d , off1[0] , off2[0] , true , false , childParent) ,
		integrator.dot(d , off1[1] , off2[1] , true , false , childParent) ,
		integrator.dot(d , off1[2] , off2[2] , true , false , childParent)
	};
	return  -Point3D< double >(dv[0] * vv[1] * vv[2], vv[0] * dv[1] * vv[2], vv[0] * vv[1] * dv[2]);
#endif // GRADIENT_DOMAIN_SOLUTION
}
template< class Real >
Point3D< double > Octree2< Real >::GetDivergence2(const typename BSplineData< 2 >::Integrator& integrator, int d, const int off1[], const int off2[], bool childParent) const
{
	double vv[] =
	{
		integrator.dot(d , off1[0] , off2[0] , false , false , childParent) ,
		integrator.dot(d , off1[1] , off2[1] , false , false , childParent) ,
		integrator.dot(d , off1[2] , off2[2] , false , false , childParent)
	};
#if GRADIENT_DOMAIN_SOLUTION
	// Take the dot-product of the vector-field with the gradient of the basis function
	double dv[] =
	{
		integrator.dot(d , off1[0] , off2[0] , true , false , childParent) ,
		integrator.dot(d , off1[1] , off2[1] , true , false , childParent) ,
		integrator.dot(d , off1[2] , off2[2] , true , false , childParent)
	};
	return  Point3D< double >(dv[0] * vv[1] * vv[2], vv[0] * dv[1] * vv[2], vv[0] * vv[1] * dv[2]);
#else // !GRADIENT_DOMAIN_SOLUTION
	// Take the dot-product of the divergence of the vector-field with the basis function
	double vd[] =
	{
		integrator.dot(d , off1[0] , off2[0] , false , true , childParent) ,
		integrator.dot(d , off1[1] , off2[1] , false , true , childParent) ,
		integrator.dot(d , off1[2] , off2[2] , false , true , childParent)
	};
	return -Point3D< double >(vd[0] * vv[1] * vv[2], vv[0] * vd[1] * vv[2], vv[0] * vv[1] * vd[2]);
#endif // GRADIENT_DOMAIN_SOLUTION
}

template< class Real >
int Octree2< Real >::GetMatrixRowSize(const typename TreeOctNode::Neighbors5& neighbors5, bool symmetric) const
{
	int count = 0;
	int nodeIndex = neighbors5.neighbors[2][2][2]->nodeData.nodeIndex;
	const TreeOctNode* const* _nodes = &neighbors5.neighbors[0][0][0];
	if (symmetric)
	{
		for (int i = 0; i < 125; i++) if (_nodes[i] && _nodes[i]->nodeData.nodeIndex >= nodeIndex) count++;
	}
	else
	{
		for (int i = 0; i < 125; i++) if (_nodes[i]) count++;
	}
	return count;
}

template< class Real >
int Octree2< Real >::SetMatrixRow(const PointInfo& pointInfo, const typename TreeOctNode::Neighbors5& neighbors5, Pointer(MatrixEntry< Real >) rowM, Pointer(MatrixEntry< Real >) rowU, Pointer(MatrixEntry< Real >) rowH, int &H_size, int offset, const typename BSplineData< 2 >::Integrator& integrator, Stencil< double, 5 > stencil, bool symmetric) const
{
	const std::vector< _PointData >& points = pointInfo.points;
	bool hasYZPoints[3], hasZPoints[3][3];
	Real diagonal = 0;
	Real splineValues[3 * 3 * 3 * 3 * 3];
	memset(splineValues, 0, sizeof(Real) * 3 * 3 * 3 * 3 * 3);

	int count = 0;
	const TreeOctNode* node = neighbors5.neighbors[2][2][2];

	bool isInterior;
	int d, off[3];
	node->depthAndOffset(d, off);
	int o = _boundaryType == 0 ? (1 << (d - 2)) : 0;
	int mn = 2 + o, mx = (1 << d) - 2 - o;
	isInterior = (off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx);

	if (_constrainValues)
	{
		int idx[3]; node->centerIndex(idx);
		for (int j = 0; j < 3; j++)
		{
			hasYZPoints[j] = false;
			for (int k = 0; k < 3; k++)
			{
				hasZPoints[j][k] = false;
				for (int l = 0; l < 3; l++)
				{
					//idx是node在x,y,z方向上的index([0, 2^d - 1])
					//_fData.baseBSplines[idx][k]就是第idx处B样条函数的的第k个小段
					const TreeOctNode* _node = neighbors5.neighbors[j + 1][k + 1][l + 1];
					if (_node && pointInfo.pointIndex(_node) != -1)
					{
						const _PointData& pData = points[pointInfo.pointIndex(_node)];
						Real* _splineValues = splineValues + 3 * 3 * (3 * (3 * j + k) + l);
						Real weight = pData.weight;
						Point3D< Real > p = pData.position;
						for (int s = 0; s < 3; s++)
						{
#if ROBERTO_TOLDO_FIX
							if (idx[0] + j - s >= 0 && idx[0] + j - s < ((2 << d) - 1)) _splineValues[3 * 0 + s] = Real(_fData.baseBSplines[idx[0] + j - s][s](p[0]));
							if (idx[1] + k - s >= 0 && idx[1] + k - s < ((2 << d) - 1)) _splineValues[3 * 1 + s] = Real(_fData.baseBSplines[idx[1] + k - s][s](p[1]));
							if (idx[2] + l - s >= 0 && idx[2] + l - s < ((2 << d) - 1)) _splineValues[3 * 2 + s] = Real(_fData.baseBSplines[idx[2] + l - s][s](p[2]));
#else // !ROBERTO_TOLDO_FIX
							_splineValues[3 * 0 + s] = Real(_fData.baseBSplines[idx[0] + j - s][s](p[0]));
							_splineValues[3 * 1 + s] = Real(_fData.baseBSplines[idx[1] + k - s][s](p[1]));
							_splineValues[3 * 2 + s] = Real(_fData.baseBSplines[idx[2] + l - s][s](p[2]));
#endif // ROBERTO_TOLDO_FIX
						}
						Real value = _splineValues[3 * 0 + j] * _splineValues[3 * 1 + k] * _splineValues[3 * 2 + l];
						Real weightedValue = value * weight;
						for (int s = 0; s < 3; s++) _splineValues[3 * 0 + s] *= weightedValue;
						diagonal += value * value * weight;
						hasYZPoints[j] = hasZPoints[j][k] = true;
					}
				}
			}
		}
	}

	Real pointValues[5][5][5];
	if (_constrainValues)
	{
		memset(pointValues, 0, sizeof(Real) * 5 * 5 * 5);
		for (int i = 0; i < 3; i++) if (hasYZPoints[i])
			for (int j = 0; j < 3; j++) if (hasZPoints[i][j])
				for (int k = 0; k < 3; k++)
				{
					const TreeOctNode* _node = neighbors5.neighbors[i + 1][j + 1][k + 1];
					if (_node && pointInfo.pointIndex(_node) != -1)
					{
						const Real* _splineValuesX = splineValues + 3 * (3 * (3 * (3 * i + j) + k) + 0) + 2;
						const Real* _splineValuesY = splineValues + 3 * (3 * (3 * (3 * i + j) + k) + 1) + 2;
						const Real* _splineValuesZ = splineValues + 3 * (3 * (3 * (3 * i + j) + k) + 2) + 2;
						for (int ii = 0; ii <= 2; ii++)
						{
							Real splineValue = _splineValuesX[-ii];
							for (int jj = 0; jj <= 2; jj++)
							{
								Real* _pointValues = pointValues[i + ii][j + jj] + k;
								Real _splineValue = splineValue * _splineValuesY[-jj];
								for (int kk = 0; kk <= 2; kk++) _pointValues[kk] += _splineValue * _splineValuesZ[-kk];
							}
						}
					}
				}
	}

	pointValues[2][2][2] = diagonal;
	int nodeIndex = neighbors5.neighbors[2][2][2]->nodeData.nodeIndex;
	if (isInterior) // General case, so try to make fast
	{
		const TreeOctNode* const* _nodes = &neighbors5.neighbors[0][0][0];
		const double* _stencil = &stencil.values[0][0][0];
		Real* _values = &pointValues[0][0][0];
		if (symmetric)
		{
			pointValues[2][2][2] /= 2;
		}
		rowM[count] = MatrixEntry< Real >(nodeIndex - offset, _stencil[5 * 5 * 2 + 5 * 2 + 2] / 2);
		rowU[count] = MatrixEntry< Real >(nodeIndex - offset, _values[5 * 5 * 2 + 5 * 2 + 2]);
		count++;
		if (symmetric)
		{
			for (int i = 0; i < 125; i++)
			{
				if (i != (5 * 5 * 2 + 5 * 2 + 2) && _nodes[i] && _nodes[i]->nodeData.nodeIndex >= nodeIndex)
				{
					rowM[count] = MatrixEntry< Real >(_nodes[i]->nodeData.nodeIndex - offset, _stencil[i]);
					rowU[count] = MatrixEntry< Real >(_nodes[i]->nodeData.nodeIndex - offset, _values[i]);
					count++;
				}
			}
		}
		else
		{
			for (int i = 0; i < 125; i++) if (i != (5 * 5 * 2 + 5 * 2 + 2) && _nodes[i])
			{
				rowM[count] = MatrixEntry< Real >(_nodes[i]->nodeData.nodeIndex - offset, _stencil[i]);
				rowU[count] = MatrixEntry< Real >(_nodes[i]->nodeData.nodeIndex - offset, _values[i]);
				count++;
			}
		}
	}
	else
	{
		int d, off[3];
		node->depthAndOffset(d, off);
		Real temp1 = Real(GetLaplacian(integrator, d, off, off, false));
		Real temp2 = pointValues[2][2][2];
		if (symmetric)
		{
			temp1 /= 2;
			temp2 /= 2;
		}
		rowM[count] = MatrixEntry< Real >(nodeIndex - offset, temp1);
		rowU[count] = MatrixEntry< Real >(nodeIndex - offset, temp2);
		count++;
		for (int x = 0; x < 5; x++) for (int y = 0; y < 5; y++) for (int z = 0; z < 5; z++)
		{
			if ((x != 2 || y != 2 || z != 2) && neighbors5.neighbors[x][y][z] && neighbors5.neighbors[x][y][z]->nodeData.nodeIndex >= 0 && (!symmetric || neighbors5.neighbors[x][y][z]->nodeData.nodeIndex >= nodeIndex))
			{
				const TreeOctNode* _node = neighbors5.neighbors[x][y][z];
				int _d, _off[3];
				_node->depthAndOffset(_d, _off);
				Real temp1 = Real(GetLaplacian(integrator, d, off, _off, false));
				Real temp2 = pointValues[x][y][z];
				if (symmetric && x == 2 && y == 2 && z == 2)
				{
					temp1 /= 2;
					temp2 /= 2;
				}
				rowM[count] = MatrixEntry< Real >(_node->nodeData.nodeIndex - offset, temp1);
				rowU[count] = MatrixEntry< Real >(_node->nodeData.nodeIndex - offset, temp2);
				count++;
			}
		}
	}
	return count;
}
template< class Real >
int Octree2< Real >::SetMatrixRow2(const PointInfo& pointInfo, const typename TreeOctNode::Neighbors5& neighbors5, Pointer(MatrixEntry< Real >) rowM, int offset, const typename BSplineData< 2 >::Integrator& integrator, Stencil< double, 5 > stencil, bool symmetric) const
{
	/*
	const std::vector< _PointData >& points = pointInfo.points;
	bool hasYZPoints[3], hasZPoints[3][3];
	Real diagonal = 0;
	Real splineValues[3 * 3 * 3 * 3 * 3];
	memset(splineValues, 0, sizeof(Real) * 3 * 3 * 3 * 3 * 3);
	*/
	int count = 0;
	const TreeOctNode* node = neighbors5.neighbors[2][2][2];

	bool isInterior;
	int d, off[3];
	node->depthAndOffset(d, off);
	int o = _boundaryType == 0 ? (1 << (d - 2)) : 0;
	int mn = 2 + o, mx = (1 << d) - 2 - o;
	isInterior = (off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx);
	/*
	if (_constrainValues)
	{
		int idx[3]; node->centerIndex(idx);
		for (int j = 0; j < 3; j++)
		{
			hasYZPoints[j] = false;
			for (int k = 0; k < 3; k++)
			{
				hasZPoints[j][k] = false;
				for (int l = 0; l < 3; l++)
				{
					//idx是node在x,y,z方向上的index([0, 2^d - 1])
					//_fData.baseBSplines[idx][k]就是第idx处B样条函数的的第k个小段
					const TreeOctNode* _node = neighbors5.neighbors[j + 1][k + 1][l + 1];
					if (_node && pointInfo.pointIndex(_node) != -1)
					{
						const _PointData& pData = points[pointInfo.pointIndex(_node)];
						Real* _splineValues = splineValues + 3 * 3 * (3 * (3 * j + k) + l);
						Real weight = pData.weight;
						Point3D< Real > p = pData.position;
						for (int s = 0; s < 3; s++)
						{
#if ROBERTO_TOLDO_FIX
							if (idx[0] + j - s >= 0 && idx[0] + j - s < ((2 << d) - 1)) _splineValues[3 * 0 + s] = Real(_fData.baseBSplines[idx[0] + j - s][s](p[0]));
							if (idx[1] + k - s >= 0 && idx[1] + k - s < ((2 << d) - 1)) _splineValues[3 * 1 + s] = Real(_fData.baseBSplines[idx[1] + k - s][s](p[1]));
							if (idx[2] + l - s >= 0 && idx[2] + l - s < ((2 << d) - 1)) _splineValues[3 * 2 + s] = Real(_fData.baseBSplines[idx[2] + l - s][s](p[2]));
#else // !ROBERTO_TOLDO_FIX
							_splineValues[3 * 0 + s] = Real(_fData.baseBSplines[idx[0] + j - s][s](p[0]));
							_splineValues[3 * 1 + s] = Real(_fData.baseBSplines[idx[1] + k - s][s](p[1]));
							_splineValues[3 * 2 + s] = Real(_fData.baseBSplines[idx[2] + l - s][s](p[2]));
#endif // ROBERTO_TOLDO_FIX
						}
						Real value = _splineValues[3 * 0 + j] * _splineValues[3 * 1 + k] * _splineValues[3 * 2 + l];
						Real weightedValue = value * weight;
						for (int s = 0; s < 3; s++) _splineValues[3 * 0 + s] *= weightedValue;
						diagonal += value * value * weight;
						hasYZPoints[j] = hasZPoints[j][k] = true;
					}
				}
			}
		}
	}
	
	Real pointValues[5][5][5];
	if (_constrainValues)
	{
		memset(pointValues, 0, sizeof(Real) * 5 * 5 * 5);
		for (int i = 0; i < 3; i++) if (hasYZPoints[i])
			for (int j = 0; j < 3; j++) if (hasZPoints[i][j])
				for (int k = 0; k < 3; k++)
				{
					const TreeOctNode* _node = neighbors5.neighbors[i + 1][j + 1][k + 1];
					if (_node && pointInfo.pointIndex(_node) != -1)
					{
						const Real* _splineValuesX = splineValues + 3 * (3 * (3 * (3 * i + j) + k) + 0) + 2;
						const Real* _splineValuesY = splineValues + 3 * (3 * (3 * (3 * i + j) + k) + 1) + 2;
						const Real* _splineValuesZ = splineValues + 3 * (3 * (3 * (3 * i + j) + k) + 2) + 2;
						for (int ii = 0; ii <= 2; ii++)
						{
							Real splineValue = _splineValuesX[-ii];
							for (int jj = 0; jj <= 2; jj++)
							{
								Real* _pointValues = pointValues[i + ii][j + jj] + k;
								Real _splineValue = splineValue * _splineValuesY[-jj];
								for (int kk = 0; kk <= 2; kk++) _pointValues[kk] += _splineValue * _splineValuesZ[-kk];
							}
						}
					}
				}
	}
	pointValues[2][2][2] = diagonal;
	*/
	int nodeIndex = neighbors5.neighbors[2][2][2]->nodeData.nodeIndex;
	if (isInterior) // General case, so try to make fast
	{
		const TreeOctNode* const* _nodes = &neighbors5.neighbors[0][0][0];
		const double* _stencil = &stencil.values[0][0][0];
		//Real* _values = &pointValues[0][0][0];
		/*
		if (symmetric)
		{
			pointValues[2][2][2] /= 2;
		}
		*/
		rowM[count] = MatrixEntry< Real >(nodeIndex - offset, _stencil[5 * 5 * 2 + 5 * 2 + 2] / 2);
		count++;
		if (symmetric)
		{
			for (int i = 0; i < 125; i++)
			{
				if (i != (5 * 5 * 2 + 5 * 2 + 2) && _nodes[i] && _nodes[i]->nodeData.nodeIndex >= nodeIndex)
				{
					rowM[count] = MatrixEntry< Real >(_nodes[i]->nodeData.nodeIndex - offset, _stencil[i]);
					count++;
				}
			}
		}
		else
		{
			for (int i = 0; i < 125; i++) if (i != (5 * 5 * 2 + 5 * 2 + 2) && _nodes[i])
			{
				rowM[count] = MatrixEntry< Real >(_nodes[i]->nodeData.nodeIndex - offset, _stencil[i]);
				count++;
			}
		}
	}
	else
	{
		int d, off[3];
		node->depthAndOffset(d, off);
		Real temp1 = Real(GetLaplacian(integrator, d, off, off, false));
		//Real temp2 = pointValues[2][2][2];
		if (symmetric)
		{
			temp1 /= 2;
			//temp2 /= 2;
		}
		rowM[count] = MatrixEntry< Real >(nodeIndex - offset, temp1);
		count++;
		for (int x = 0; x < 5; x++) for (int y = 0; y < 5; y++) for (int z = 0; z < 5; z++)
		{
			if ((x != 2 || y != 2 || z != 2) && neighbors5.neighbors[x][y][z] && neighbors5.neighbors[x][y][z]->nodeData.nodeIndex >= 0 && (!symmetric || neighbors5.neighbors[x][y][z]->nodeData.nodeIndex >= nodeIndex))
			{
				const TreeOctNode* _node = neighbors5.neighbors[x][y][z];
				int _d, _off[3];
				_node->depthAndOffset(_d, _off);
				Real temp1 = Real(GetLaplacian(integrator, d, off, _off, false));
				//Real temp2 = pointValues[x][y][z];
				if (symmetric && x == 2 && y == 2 && z == 2)
				{
					temp1 /= 2;
					//temp2 /= 2;
				}
				rowM[count] = MatrixEntry< Real >(_node->nodeData.nodeIndex - offset, temp1);
				count++;
			}
		}
	}
	return count;
}
template< class Real >
int Octree2< Real >::SetMatrixRowSPR(const PointInfo& pointInfo, const typename TreeOctNode::Neighbors5& neighbors5, Pointer(MatrixEntry< Real >) row, int offset, const typename BSplineData< 2 >::Integrator& integrator, const Stencil< double, 5 >& stencil, bool symmetric) const
{
	//SetMatrixRow求的是矩阵元素的第i行，neighbors5也是编号i的结点对应的neighbour,222就是自己
	/*
	const std::vector< _PointData >& points = pointInfo.points;
	bool hasYZPoints[3], hasZPoints[3][3];
	Real diagonal = 0;
	Real splineValues[3 * 3 * 3 * 3 * 3];
	memset(splineValues, 0, sizeof(Real) * 3 * 3 * 3 * 3 * 3);
	*/
	Real diagonal = 0;
	int count = 0;
	//i号结点自己
	const TreeOctNode* node = neighbors5.neighbors[2][2][2];

	bool isInterior;
	int d, off[3];
	node->depthAndOffset(d, off);

	int o = _boundaryType == 0 ? (1 << (d - 2)) : 0;
	int mn = 2 + o, mx = (1 << d) - 2 - o;
	isInterior = (off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx);

	/*
	if (_constrainValues)
	{
		//idx为node i在x,y,z方向的index
		int idx[3]; node->centerIndex(idx);
		for (int j = 0; j < 3; j++)
		{
			hasYZPoints[j] = false;
			for (int k = 0; k < 3; k++)
			{
				hasZPoints[j][k] = false;
				for (int l = 0; l < 3; l++)
				{
					//j k l各三维，找x,y,z三个方向的neighbour
					const TreeOctNode* _node = neighbors5.neighbors[j + 1][k + 1][l + 1];
					if (_node && pointInfo.pointIndex(_node) != -1)
					{
						const _PointData& pData = points[pointInfo.pointIndex(_node)];
						Real* _splineValues = splineValues + 3 * 3 * (3 * (3 * j + k) + l);
						Real weight = pData.weight;
						Point3D< Real > p = pData.position;
						//s代表B样条的某一段
						//idx[0]+j-s处的基函数的从左向右第s段会对idx[0]+j有贡献
						for (int s = 0; s < 3; s++)
						{
#if ROBERTO_TOLDO_FIX
							if (idx[0] + j - s >= 0 && idx[0] + j - s < ((2 << d) - 1)) _splineValues[3 * 0 + s] = Real(_fData.baseBSplines[idx[0] + j - s][s](p[0]));
							if (idx[1] + k - s >= 0 && idx[1] + k - s < ((2 << d) - 1)) _splineValues[3 * 1 + s] = Real(_fData.baseBSplines[idx[1] + k - s][s](p[1]));
							if (idx[2] + l - s >= 0 && idx[2] + l - s < ((2 << d) - 1)) _splineValues[3 * 2 + s] = Real(_fData.baseBSplines[idx[2] + l - s][s](p[2]));
#else // !ROBERTO_TOLDO_FIX
							_splineValues[3 * 0 + s] = Real(_fData.baseBSplines[idx[0] + j - s][s](p[0]));
							_splineValues[3 * 1 + s] = Real(_fData.baseBSplines[idx[1] + k - s][s](p[1]));
							_splineValues[3 * 2 + s] = Real(_fData.baseBSplines[idx[2] + l - s][s](p[2]));
#endif // ROBERTO_TOLDO_FIX
						}
						//value的值是j k l处的B样条基函数值在node(222)处的贡献
						Real value = _splineValues[3 * 0 + j] * _splineValues[3 * 1 + k] * _splineValues[3 * 2 + l];
						Real weightedValue = value * weight;
						//Bi * Bj
						for (int s = 0; s < 3; s++) _splineValues[3 * 0 + s] *= weightedValue;
						diagonal += value * value * weight;
						hasYZPoints[j] = hasZPoints[j][k] = true;
					}
				}
			}
		}
	}
	*/
	//pointValues正式存放sum(2^d <Bi, Bj>)
	Real pointValues[5][5][5];
	/*
	if (_constrainValues)
	{
		memset(pointValues, 0, sizeof(Real) * 5 * 5 * 5);
		for (int i = 0; i < 3; i++) if (hasYZPoints[i])
			for (int j = 0; j < 3; j++) if (hasZPoints[i][j])
				for (int k = 0; k < 3; k++)
				{
					const TreeOctNode* _node = neighbors5.neighbors[i + 1][j + 1][k + 1];
					if (_node && pointInfo.pointIndex(_node) != -1)
					{
						const Real* _splineValuesX = splineValues + 3 * (3 * (3 * (3 * i + j) + k) + 0) + 2;
						const Real* _splineValuesY = splineValues + 3 * (3 * (3 * (3 * i + j) + k) + 1) + 2;
						const Real* _splineValuesZ = splineValues + 3 * (3 * (3 * (3 * i + j) + k) + 2) + 2;
						for (int ii = 0; ii <= 2; ii++)
						{
							Real splineValue = _splineValuesX[-ii];
							for (int jj = 0; jj <= 2; jj++)
							{
								Real* _pointValues = pointValues[i + ii][j + jj] + k;
								Real _splineValue = splineValue * _splineValuesY[-jj];
								for (int kk = 0; kk <= 2; kk++) _pointValues[kk] += _splineValue * _splineValuesZ[-kk];
							}
						}
					}
				}
	}
	*/
	pointValues[2][2][2] = diagonal;
	int nodeIndex = neighbors5.neighbors[2][2][2]->nodeData.nodeIndex;
	if (isInterior) // General case, so try to make fast
	{
		const TreeOctNode* const* _nodes = &neighbors5.neighbors[0][0][0];
		const double* _stencil = &stencil.values[0][0][0];
		Real* _values = &pointValues[0][0][0];
		if (_constrainValues) for (int i = 0; i < 125; i++) _values[i] = Real(_stencil[i] + _values[i]);
		else                   for (int i = 0; i < 125; i++) _values[i] = Real(_stencil[i]);
		if (symmetric) pointValues[2][2][2] /= 2;
		row[count++] = MatrixEntry< Real >(nodeIndex - offset, _values[5 * 5 * 2 + 5 * 2 + 2]);
		if (symmetric)
		{
			for (int i = 0; i < 125; i++) if (i != (5 * 5 * 2 + 5 * 2 + 2) && _nodes[i] && _nodes[i]->nodeData.nodeIndex >= nodeIndex)
				row[count++] = MatrixEntry< Real >(_nodes[i]->nodeData.nodeIndex - offset, _values[i]);
		}
		else
		{
			for (int i = 0; i < 125; i++) if (i != (5 * 5 * 2 + 5 * 2 + 2) && _nodes[i])
				row[count++] = MatrixEntry< Real >(_nodes[i]->nodeData.nodeIndex - offset, _values[i]);
		}
	}
	else
	{
		int d, off[3];
		node->depthAndOffset(d, off);
		Real temp = Real(GetLaplacian(integrator, d, off, off, false));
		if (_constrainValues) temp += pointValues[2][2][2];
		if (symmetric) temp /= 2;
		row[count++] = MatrixEntry< Real >(nodeIndex - offset, temp);
		for (int x = 0; x < 5; x++) for (int y = 0; y < 5; y++) for (int z = 0; z < 5; z++)
			if ((x != 2 || y != 2 || z != 2) && neighbors5.neighbors[x][y][z] && neighbors5.neighbors[x][y][z]->nodeData.nodeIndex >= 0 && (!symmetric || neighbors5.neighbors[x][y][z]->nodeData.nodeIndex >= nodeIndex))
			{
				const TreeOctNode* _node = neighbors5.neighbors[x][y][z];
				int _d, _off[3];
				_node->depthAndOffset(_d, _off);
				Real temp = Real(GetLaplacian(integrator, d, off, _off, false));
				if (_constrainValues) temp += pointValues[x][y][z];
				if (symmetric && x == 2 && y == 2 && z == 2) temp /= 2;
				row[count++] = MatrixEntry< Real >(_node->nodeData.nodeIndex - offset, temp);
			}
	}
	return count;
}
// if( scatter ) normals come from the center node
// else          normals come from the neighbors
template< class Real >
void Octree2< Real >::SetDivergenceStencil(int depth, const typename BSplineData< 2 >::Integrator& integrator, Stencil< Point3D< double >, 5 >& stencil, bool scatter) const
{
	if (depth < 2) return;
	int center = 1 << (depth - 1);
	int offset[] = { center , center , center };
	for (int x = 0; x < 5; x++) for (int y = 0; y < 5; y++) for (int z = 0; z < 5; z++)
	{
		int _offset[] = { x + center - 2 , y + center - 2 , z + center - 2 };
		if (scatter) stencil.values[x][y][z] = GetDivergence1(integrator, depth, offset, _offset, false);
		else          stencil.values[x][y][z] = GetDivergence2(integrator, depth, offset, _offset, false);
	}
}
template< class Real >
void Octree2< Real >::SetDivergenceStencils(int depth, const typename BSplineData< 2 >::Integrator& integrator, Stencil< Point3D< double >, 5 > stencils[2][2][2], bool scatter) const
{
	if (depth < 2) return;
	int center = 1 << (depth - 1);
	for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) for (int k = 0; k < 2; k++)
	{
		int offset[] = { center + i , center + j , center + k };
		for (int x = 0; x < 5; x++) for (int y = 0; y < 5; y++) for (int z = 0; z < 5; z++)
		{
			int _offset[] = { x - 2 + center / 2 , y - 2 + center / 2 , z - 2 + center / 2 };
			if (scatter) stencils[i][j][k].values[x][y][z] = GetDivergence1(integrator, depth, offset, _offset, true);
			else          stencils[i][j][k].values[x][y][z] = GetDivergence2(integrator, depth, offset, _offset, true);
		}
	}
}
template< class Real >
void Octree2< Real >::SetLaplacianStencil(int depth, const typename BSplineData< 2 >::Integrator& integrator, Stencil< double, 5 >& stencil) const
{
	if (depth < 2) return;
	int center = 1 << (depth - 1);
	int offset[] = { center , center , center };
	for (int x = -2; x <= 2; x++) for (int y = -2; y <= 2; y++) for (int z = -2; z <= 2; z++)
	{
		int _offset[] = { x + center , y + center , z + center };
		stencil.values[x + 2][y + 2][z + 2] = GetLaplacian(integrator, depth, offset, _offset, false);
	}
}
template< class Real >
void Octree2< Real >::SetLaplacianStencils(int depth, const typename BSplineData< 2 >::Integrator& integrator, Stencil< double, 5 > stencils[2][2][2]) const
{
	if (depth < 2) return;
	int center = 1 << (depth - 1);
	for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) for (int k = 0; k < 2; k++)
	{
		int offset[] = { center + i , center + j , center + k };
		for (int x = -2; x <= 2; x++) for (int y = -2; y <= 2; y++) for (int z = -2; z <= 2; z++)
		{
			int _offset[] = { x + center / 2 , y + center / 2 , z + center / 2 };
			stencils[i][j][k].values[x + 2][y + 2][z + 2] = GetLaplacian(integrator, depth, offset, _offset, true);
		}
	}
}
template< class Real >
void Octree2< Real >::SetCenterEvaluationStencil(const typename BSplineData< 2 >::template CenterEvaluator< 1 >& evaluator, int depth, Stencil< double, 3 >& stencil) const
{
	if (depth < 2) return;
	int center = 1 << (depth - 1);
	for (int x = 0; x < 3; x++) for (int y = 0; y < 3; y++) for (int z = 0; z < 3; z++)
	{
		int off[] = { center + x - 1 , center + y - 1 , center + z - 1 };
		stencil.values[x][y][z] = Real(evaluator.value(depth, center, off[0], false, false) * evaluator.value(depth, center, off[1], false, false) * evaluator.value(depth, center, off[2], false, false));
	}
}
template< class Real >
void Octree2< Real >::SetCenterEvaluationStencils(const typename BSplineData< 2 >::template CenterEvaluator< 1 >& evaluator, int depth, Stencil< double, 3 > stencils[8]) const
{
	if (depth < 3) return;
	int center = 1 << (depth - 1);
	for (int cx = 0; cx < 2; cx++) for (int cy = 0; cy < 2; cy++) for (int cz = 0; cz < 2; cz++)
	{
		int idx[] = { center + cx , center + cy , center + cz };
		for (int x = 0; x < 3; x++) for (int y = 0; y < 3; y++) for (int z = 0; z < 3; z++)
		{
			int off[] = { center / 2 + x - 1 , center / 2 + y - 1 , center / 2 + z - 1 };
			stencils[Cube::CornerIndex(cx, cy, cz)].values[x][y][z] = Real(evaluator.value(depth, idx[0], off[0], false, true) * evaluator.value(depth, idx[1], off[1], false, true) * evaluator.value(depth, idx[2], off[2], false, true));
		}
	}
}
template< class Real >
void Octree2< Real >::SetCornerEvaluationStencil(const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator, int depth, Stencil< double, 3 > stencil[8]) const
{
	if (depth < 2) return;
	int center = 1 << (depth - 1);
	for (int cx = 0; cx < 2; cx++) for (int cy = 0; cy < 2; cy++) for (int cz = 0; cz < 2; cz++)
	{
		int c = Cube::CornerIndex(cx, cy, cz);
		for (int x = 0; x < 3; x++) for (int y = 0; y < 3; y++) for (int z = 0; z < 3; z++)
		{
			int off[] = { center + x - 1 , center + y - 1 , center + z - 1 };
			stencil[c].values[x][y][z] = evaluator.value(depth, center, cx, off[0], false, false) * evaluator.value(depth, center, cy, off[1], false, false) * evaluator.value(depth, center, cz, off[2], false, false);
		}
	}
}
template< class Real >
void Octree2< Real >::SetCornerEvaluationStencils(const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator, int depth, Stencil< double, 3 > stencils[8][8]) const
{
	if (depth < 3) return;
	int center = 1 << (depth - 1);
	for (int cx = 0; cx < 2; cx++) for (int cy = 0; cy < 2; cy++) for (int cz = 0; cz < 2; cz++)
	{
		int c = Cube::CornerIndex(cx, cy, cz);
		for (int _cx = 0; _cx < 2; _cx++) for (int _cy = 0; _cy < 2; _cy++) for (int _cz = 0; _cz < 2; _cz++)
		{
			int _c = Cube::CornerIndex(_cx, _cy, _cz);
			int idx[] = { center + _cx , center + _cy , center + _cz };
			for (int x = 0; x < 3; x++) for (int y = 0; y < 3; y++) for (int z = 0; z < 3; z++)
			{
				int off[] = { center / 2 + x - 1 , center / 2 + y - 1 , center / 2 + z - 1 };
				stencils[c][_c].values[x][y][z] = evaluator.value(depth, idx[0], cx, off[0], false, true) * evaluator.value(depth, idx[1], cy, off[1], false, true) * evaluator.value(depth, idx[2], cz, off[2], false, true);
			}
		}
	}
}
template< class Real >
void Octree2< Real >::SetCornerNormalEvaluationStencil(const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator, int depth, Stencil< Point3D< double >, 3 > stencil[8]) const
{
	if (depth < 2) return;
	int center = 1 << (depth - 1);
	for (int cx = 0; cx < 2; cx++) for (int cy = 0; cy < 2; cy++) for (int cz = 0; cz < 2; cz++)
	{
		int c = Cube::CornerIndex(cx, cy, cz);
		for (int x = 0; x < 3; x++) for (int y = 0; y < 3; y++) for (int z = 0; z < 3; z++)
		{
			int off[] = { center + x - 1 , center + y - 1 , center + z - 1 };
			double v[] = { evaluator.value(depth , center , cx , off[0] , false , false) , evaluator.value(depth , center , cy , off[1] , false , false) , evaluator.value(depth , center , cz , off[2] , false , false) };
			double dv[] = { evaluator.value(depth , center , cx , off[0] , true  , false) , evaluator.value(depth , center , cy , off[1] , true  , false) , evaluator.value(depth , center , cz , off[2] , true  , false) };
			stencil[c].values[x][y][z] = Point3D< double >(dv[0] * v[1] * v[2], v[0] * dv[1] * v[2], v[0] * v[1] * dv[2]);
		}
	}
}
template< class Real >
void Octree2< Real >::SetCornerNormalEvaluationStencils(const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator, int depth, Stencil< Point3D< double >, 3 > stencils[8][8]) const
{
	if (depth < 3) return;
	int center = 1 << (depth - 1);
	for (int cx = 0; cx < 2; cx++) for (int cy = 0; cy < 2; cy++) for (int cz = 0; cz < 2; cz++)
	{
		int c = Cube::CornerIndex(cx, cy, cz);	// Which corner of the finer cube
		for (int _cx = 0; _cx < 2; _cx++) for (int _cy = 0; _cy < 2; _cy++) for (int _cz = 0; _cz < 2; _cz++)
		{
			int _c = Cube::CornerIndex(_cx, _cy, _cz);	// Which child node
			int idx[] = { center + _cx , center + _cy , center + _cz };
			for (int x = 0; x < 3; x++) for (int y = 0; y < 3; y++) for (int z = 0; z < 3; z++)
			{
				int off[] = { center / 2 + x - 1 , center / 2 + y - 1 , center / 2 + z - 1 };
				double v[] = { evaluator.value(depth , idx[0] , cx , off[0] , false , true) , evaluator.value(depth , idx[1] , cy , off[1] , false , true) , evaluator.value(depth , idx[2] , cz , off[2] , false , true) };
				double dv[] = { evaluator.value(depth , idx[0] , cx , off[0] , true  , true) , evaluator.value(depth , idx[1] , cy , off[1] , true  , true) , evaluator.value(depth , idx[2] , cz , off[2] , true  , true) };
				stencils[c][_c].values[x][y][z] = Point3D< double >(dv[0] * v[1] * v[2], v[0] * dv[1] * v[2], v[0] * v[1] * dv[2]);
			}
		}
	}
}
template< class Real >
void Octree2< Real >::SetCornerNormalEvaluationStencil(const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator, int depth, Stencil< Point3D< double >, 5 > stencil[8]) const
{
	if (depth < 2) return;
	int center = 1 << (depth - 1);
	for (int cx = 0; cx < 2; cx++) for (int cy = 0; cy < 2; cy++) for (int cz = 0; cz < 2; cz++)
	{
		int c = Cube::CornerIndex(cx, cy, cz);
		for (int x = 0; x < 5; x++) for (int y = 0; y < 5; y++) for (int z = 0; z < 5; z++)
		{
			int off[] = { center + x - 2 , center + y - 2 , center + z - 2 };
			double v[] = { evaluator.value(depth , center , cx , off[0] , false , false) , evaluator.value(depth , center , cy , off[1] , false , false) , evaluator.value(depth , center , cz , off[2] , false , false) };
			double dv[] = { evaluator.value(depth , center , cx , off[0] , true  , false) , evaluator.value(depth , center , cy , off[1] , true  , false) , evaluator.value(depth , center , cz , off[2] , true  , false) };
			stencil[c].values[x][y][z] = Point3D< double >(dv[0] * v[1] * v[2], v[0] * dv[1] * v[2], v[0] * v[1] * dv[2]);
		}
	}
}
template< class Real >
void Octree2< Real >::SetCornerNormalEvaluationStencils(const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator, int depth, Stencil< Point3D< double >, 5 > stencils[8][8]) const
{
	if (depth < 3) return;
	int center = 1 << (depth - 1);
	for (int cx = 0; cx < 2; cx++) for (int cy = 0; cy < 2; cy++) for (int cz = 0; cz < 2; cz++)
	{
		int c = Cube::CornerIndex(cx, cy, cz);	// Which corner of the finer cube
		for (int _cx = 0; _cx < 2; _cx++) for (int _cy = 0; _cy < 2; _cy++) for (int _cz = 0; _cz < 2; _cz++)
		{
			int _c = Cube::CornerIndex(_cx, _cy, _cz);	// Which child node
			int idx[] = { center + _cx , center + _cy , center + _cz };
			for (int x = 0; x < 5; x++) for (int y = 0; y < 5; y++) for (int z = 0; z < 5; z++)
			{
				int off[] = { center / 2 + x - 2 , center / 2 + y - 2 , center / 2 + z - 2 };
				double v[] = { evaluator.value(depth , idx[0] , cx , off[0] , false , true) , evaluator.value(depth , idx[1] , cy , off[1] , false , true) , evaluator.value(depth , idx[2] , cz , off[2] , false , true) };
				double dv[] = { evaluator.value(depth , idx[0] , cx , off[0] , true  , true) , evaluator.value(depth , idx[1] , cy , off[1] , true  , true) , evaluator.value(depth , idx[2] , cz , off[2] , true  , true) };
				stencils[c][_c].values[x][y][z] = Point3D< double >(dv[0] * v[1] * v[2], v[0] * dv[1] * v[2], v[0] * v[1] * dv[2]);
			}
		}
	}
}

template< class Real >
void Octree2< Real >::UpdateCoarserSupportBounds(const TreeOctNode* node, int& startX, int& endX, int& startY, int& endY, int& startZ, int& endZ)
{
	if (node->parent)
	{
		int x, y, z, c = int(node - node->parent->children);
		Cube::FactorCornerIndex(c, x, y, z);
		if (x == 0) endX = 4;
		else     startX = 1;
		if (y == 0) endY = 4;
		else     startY = 1;
		if (z == 0) endZ = 4;
		else     startZ = 1;
	}
}
// Given the solution @( depth ) add to the met constraints @( depth-1 )
template< class Real >
void Octree2< Real >::UpdateConstraintsFromFiner(const typename BSplineData< 2 >::Integrator& integrator, int depth, const SortedTreeNodes& sNodes, ConstPointer(Real) fineSolution, Pointer(Real) coarseConstraints) const
{
	if (depth <= _minDepth) return;
	Stencil< double, 5 > stencils[2][2][2];
	// Get the stencil describing the Laplacian relating coefficients @(depth) with coefficients @(depth-1)
	SetLaplacianStencils(depth, integrator, stencils);
	size_t start = sNodes.nodeCount[depth], end = sNodes.nodeCount[depth + 1], range = end - start;
	int lStart = sNodes.nodeCount[depth - 1];
	memset(coarseConstraints, 0, sizeof(Real) * (sNodes.nodeCount[depth] - sNodes.nodeCount[depth - 1]));

	// Iterate over the nodes @( depth )
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
	for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(depth - 1);
#pragma omp parallel for num_threads( threads )
	for (int i = sNodes.nodeCount[depth]; i < sNodes.nodeCount[depth + 1]; i++)
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
		TreeOctNode* node = sNodes.treeNodes[i];

		bool insetSupported = _boundaryType != 0 || _IsInsetSupported(node);

		// Offset the coarser constraints using the solution from the current resolutions.
		int x, y, z, c;
		c = int(node - node->parent->children);
		Cube::FactorCornerIndex(c, x, y, z);
		if (insetSupported)
		{
			typename TreeOctNode::Neighbors5 pNeighbors5;
			neighborKey.getNeighbors(node->parent, pNeighbors5);
			const Stencil< double, 5 >& lapStencil = stencils[x][y][z];

			Pointer(Real) __coarseConstraints = coarseConstraints - lStart;
			bool isInterior;
			int d, off[3];
			{
				node->depthAndOffset(d, off);
				int o = _boundaryType == 0 ? (1 << (d - 2)) : 0;
				int mn = 4 + o, mx = (1 << d) - 4 - o;
				isInterior = (off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx);
			}
			// Offset the constraints using the solution from finer resolutions.
			int startX = 0, endX = 5, startY = 0, endY = 5, startZ = 0, endZ = 5;
			UpdateCoarserSupportBounds(node, startX, endX, startY, endY, startZ, endZ);

			Real solution = fineSolution[node->nodeData.nodeIndex - start];
			for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
				if (pNeighbors5.neighbors[x][y][z] && pNeighbors5.neighbors[x][y][z]->nodeData.nodeIndex >= 0)
				{
					const TreeOctNode* _node = pNeighbors5.neighbors[x][y][z];
					if (isInterior)
#pragma omp atomic
						__coarseConstraints[_node->nodeData.nodeIndex] += Real(lapStencil.values[x][y][z] * solution);
					else
					{
						int _d, _off[3];
						_node->depthAndOffset(_d, _off);
#pragma omp atomic
						__coarseConstraints[_node->nodeData.nodeIndex] += Real(GetLaplacian(integrator, d, off, _off, true) * solution);
					}
				}
		}
	}
}

template< class Real >
void Octree2< Real >::UpdateConstraintsFromCoarser(const PointInfo& pointInfo, const typename TreeOctNode::Neighbors5& neighbors5, const typename TreeOctNode::Neighbors5& pNeighbors5, TreeOctNode* node, Vector<Real> &C, ConstPointer(Real) metSolution, const typename BSplineData< 2 >::Integrator& integrator, const Stencil< double, 5 >& lapStencil) const
{
	const std::vector< _PointData >& points = pointInfo.points;
	if (node->depth() <= _minDepth) return;
	bool isInterior;
	int d, off[3];
	{
		node->depthAndOffset(d, off);
		int o = _boundaryType == 0 ? (1 << (d - 2)) : 0;
		int mn = 4 + o, mx = (1 << d) - 4 - o;
		isInterior = (off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx);
	}
	Real constraint = Real(0);
	// Offset the constraints using the solution from lower resolutions.
	int startX = 0, endX = 5, startY = 0, endY = 5, startZ = 0, endZ = 5;
	UpdateCoarserSupportBounds(node, startX, endX, startY, endY, startZ, endZ);

	for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
	{
		if (pNeighbors5.neighbors[x][y][z] && pNeighbors5.neighbors[x][y][z]->nodeData.nodeIndex >= 0)
		{
			const TreeOctNode* _node = pNeighbors5.neighbors[x][y][z];
			Real _solution = metSolution[_node->nodeData.nodeIndex];
			{
				if (isInterior) C[node->nodeData.nodeIndex] += Real(lapStencil.values[x][y][z] * _solution);
				else
				{
					int _d, _off[3];
					_node->depthAndOffset(_d, _off);
					C[node->nodeData.nodeIndex] += Real(GetLaplacian(integrator, d, off, _off, true) * _solution);
				}
			}
		}
	}
	/*
	if (_constrainValues)
	{
		double constraint = 0;
		int idx[3];
		node->centerIndex(idx);
		// Evaluate the current node's basis function at adjacent points
		for (int x = 1; x < 4; x++) for (int y = 1; y < 4; y++) for (int z = 1; z < 4; z++)
			if (neighbors5.neighbors[x][y][z] && pointInfo.pointIndex(neighbors5.neighbors[x][y][z]) != -1)
			{
				const _PointData& pData = points[pointInfo.pointIndex(neighbors5.neighbors[x][y][z])];
				Real weightedPointValue = pData.weightedCoarserValue;
				Point3D< Real > p = pData.position;
				constraint +=
					_fData.baseBSplines[idx[0]][x - 1](p[0]) *
					_fData.baseBSplines[idx[1]][y - 1](p[1]) *
					_fData.baseBSplines[idx[2]][z - 1](p[2]) *
					weightedPointValue;
			}
		C[node->nodeData.nodeIndex] += Real(constraint);
	}
	*/
}

template< class Real >
void Octree2< Real >::UpdateConstraintsFromFiner(const PointInfo& pointInfo, const typename TreeOctNode::Neighbors5& neighbors5, const typename TreeOctNode::Neighbors5& pNeighbors5, TreeOctNode* node, Vector<Real>& C, ConstPointer(Real) metSolution, const typename BSplineData< 2 >::Integrator& integrator, const Stencil< double, 5 >& lapStencil) const
{
	const std::vector< _PointData >& points = pointInfo.points;
	if (node->depth() <= _minDepth) return;
	bool isInterior;
	int d, off[3];
	{
		node->depthAndOffset(d, off);
		int o = _boundaryType == 0 ? (1 << (d - 2)) : 0;
		int mn = 4 + o, mx = (1 << d) - 4 - o;
		isInterior = (off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx);
	}
	Real constraint = Real(0);
	// Offset the constraints using the solution from lower resolutions.
	int startX = 0, endX = 5, startY = 0, endY = 5, startZ = 0, endZ = 5;
	UpdateCoarserSupportBounds(node, startX, endX, startY, endY, startZ, endZ);

	for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
	{
		if (pNeighbors5.neighbors[x][y][z] && pNeighbors5.neighbors[x][y][z]->nodeData.nodeIndex >= 0)
		{
			const TreeOctNode* _node = pNeighbors5.neighbors[x][y][z];
			Real _solution = metSolution[node->nodeData.nodeIndex];
			{
				if (isInterior) C[_node->nodeData.nodeIndex] += Real(lapStencil.values[x][y][z] * _solution);
				else
				{
					int _d, _off[3];
					_node->depthAndOffset(_d, _off);
					C[_node->nodeData.nodeIndex] += Real(GetLaplacian(integrator, d, off, _off, true) * _solution);
				}
			}
		}
	}
}


template< class Real >
void Octree2< Real >::UpdateConstraintsFromCoarserSPR(const PointInfo& pointInfo, const typename TreeOctNode::Neighbors5& neighbors5, const typename TreeOctNode::Neighbors5& pNeighbors5, TreeOctNode* node, Pointer(Real) constraints, ConstPointer(Real) metSolution, const typename BSplineData< 2 >::Integrator& integrator, const Stencil< double, 5 >& lapStencil) const
{
	const std::vector< _PointData >& points = pointInfo.points;
	if (node->depth() <= _minDepth) return;
	bool isInterior;
	int d, off[3];
	{
		node->depthAndOffset(d, off);
		int o = _boundaryType == 0 ? (1 << (d - 2)) : 0;
		int mn = 4 + o, mx = (1 << d) - 4 - o;
		isInterior = (off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx);
	}
	Real constraint = Real(0);
	// Offset the constraints using the solution from lower resolutions.
	int startX = 0, endX = 5, startY = 0, endY = 5, startZ = 0, endZ = 5;
	UpdateCoarserSupportBounds(node, startX, endX, startY, endY, startZ, endZ);

	for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
		if (pNeighbors5.neighbors[x][y][z] && pNeighbors5.neighbors[x][y][z]->nodeData.nodeIndex >= 0)
		{
			const TreeOctNode* _node = pNeighbors5.neighbors[x][y][z];
			Real _solution = metSolution[_node->nodeData.nodeIndex];
			{
				if (isInterior)
				{
					//A^{d,d-1}x^{d-1}
					constraints[node->nodeData.nodeIndex] -= Real(lapStencil.values[x][y][z] * _solution);
				}
				else
				{
					int _d, _off[3];
					_node->depthAndOffset(_d, _off);
					constraints[node->nodeData.nodeIndex] -= Real(GetLaplacian(integrator, d, off, _off, true) * _solution);
				}
			}
		}
	/*
	if (_constrainValues)
	{
		double constraint = 0;
		int idx[3];
		node->centerIndex(idx);
		// Evaluate the current node's basis function at adjacent points
		for (int x = 1; x < 4; x++) for (int y = 1; y < 4; y++) for (int z = 1; z < 4; z++)
			if (neighbors5.neighbors[x][y][z] && pointInfo.pointIndex(neighbors5.neighbors[x][y][z]) != -1)
			{
				const _PointData& pData = points[pointInfo.pointIndex(neighbors5.neighbors[x][y][z])];
				Real weightedPointValue = pData.weightedCoarserValue;
				Point3D< Real > p = pData.position;
				constraint +=
					_fData.baseBSplines[idx[0]][x - 1](p[0]) *
					_fData.baseBSplines[idx[1]][y - 1](p[1]) *
					_fData.baseBSplines[idx[2]][z - 1](p[2]) *
					weightedPointValue;
			}
		constraints[node->nodeData.nodeIndex] -= Real(constraint);
	}
	*/
}

template< class Real >
template< class C >
void Octree2< Real >::DownSample(int depth, const SortedTreeNodes& sNodes, ConstPointer(C) fineConstraints, Pointer(C) coarseConstraints) const
{
	if (depth == 0) return;
	double cornerValue;
	if (_boundaryType == -1) cornerValue = 0.50;
	else if (_boundaryType == 1) cornerValue = 1.00;
	else                         cornerValue = 0.75;
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
	for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(depth);
#pragma omp parallel for num_threads( threads )
	for (int i = sNodes.nodeCount[depth]; i < sNodes.nodeCount[depth + 1]; i++)
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
		int d, off[3];
		UpSampleData usData[3];
		sNodes.treeNodes[i]->depthAndOffset(d, off);
		for (int dd = 0; dd < 3; dd++)
		{
			if (off[dd] == 0) usData[dd] = UpSampleData(1, cornerValue, 0.00);
			else if (off[dd] + 1 == (1 << depth)) usData[dd] = UpSampleData(0, 0.00, cornerValue);
			else if (off[dd] % 2) usData[dd] = UpSampleData(1, 0.75, 0.25);
			else                             usData[dd] = UpSampleData(0, 0.25, 0.75);
		}
		typename TreeOctNode::Neighbors3& neighbors = neighborKey.getNeighbors(sNodes.treeNodes[i]->parent);
		C c = fineConstraints[i - sNodes.nodeCount[depth]];
		for (int ii = 0; ii < 2; ii++)
		{
			int _ii = ii + usData[0].start;
			C cx = C(c * usData[0].v[ii]);
			for (int jj = 0; jj < 2; jj++)
			{
				int _jj = jj + usData[1].start;
				C cxy = C(cx * usData[1].v[jj]);
				for (int kk = 0; kk < 2; kk++)
				{
					int _kk = kk + usData[2].start;
					TreeOctNode* pNode = neighbors.neighbors[_ii][_jj][_kk];
					if (pNode)
#pragma omp atomic
						coarseConstraints[pNode->nodeData.nodeIndex - sNodes.nodeCount[depth - 1]] += C(cxy * usData[2].v[kk]);
				}
			}
		}
	}
}


template< class Real >
void Octree2< Real >::DownSample2(int depth, const SortedTreeNodes& sNodes, int fineConstraints, int coarseConstraints, std::vector<std::vector<MatrixEntry <Real>>> &_Q) const
{
	if (depth == 0) return;
	double cornerValue;
	if (_boundaryType == -1) cornerValue = 0.50;
	else if (_boundaryType == 1) cornerValue = 1.00;
	else                         cornerValue = 0.75;
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
	for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(depth);
#pragma omp parallel for num_threads( threads )
	for (int i = sNodes.nodeCount[depth]; i < sNodes.nodeCount[depth + 1]; i++)
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
		int d, off[3];
		//B样条的coarse函数写成fine函数的线性组合
		UpSampleData usData[3];
		sNodes.treeNodes[i]->depthAndOffset(d, off);
		for (int dd = 0; dd < 3; dd++)
		{
			if (off[dd] == 0) usData[dd] = UpSampleData(1, cornerValue, 0.00);
			else if (off[dd] + 1 == (1 << depth)) usData[dd] = UpSampleData(0, 0.00, cornerValue);
			else if (off[dd] % 2) usData[dd] = UpSampleData(1, 0.75, 0.25);
			else                             usData[dd] = UpSampleData(0, 0.25, 0.75);
		}
		typename TreeOctNode::Neighbors3& neighbors = neighborKey.getNeighbors(sNodes.treeNodes[i]->parent);
		int fine_index = fineConstraints + i - sNodes.nodeCount[depth];
		for (int ii = 0; ii < 2; ii++)
		{
			int _ii = ii + usData[0].start;
			Real cx = Real(usData[0].v[ii]);
			for (int jj = 0; jj < 2; jj++)
			{
				int _jj = jj + usData[1].start;
				Real cxy = Real(cx * usData[1].v[jj]);
				for (int kk = 0; kk < 2; kk++)
				{
					int _kk = kk + usData[2].start;
					TreeOctNode* pNode = neighbors.neighbors[_ii][_jj][_kk];
					if (pNode)
					{
#pragma omp critical
						{
							int coarse_index = coarseConstraints + pNode->nodeData.nodeIndex - sNodes.nodeCount[depth - 1];
							AddSparseVector(_Q[coarse_index], _Q[fine_index], Real(cxy * usData[2].v[kk]));
						}
						//coarseConstraints[pNode->nodeData.nodeIndex - sNodes.nodeCount[depth - 1]] += C(cxy * usData[2].v[kk]);
					}
				}
			}
		}
	}
}

template< class Real >
void Octree2< Real >::DownSample3(int depth, const SortedTreeNodes& sNodes, int fineConstraints, int coarseConstraints, std::vector<std::vector<MatrixEntry <Real>>>& _Q, std::vector < std::unordered_map<int, int>> &_Q_index, Real thre) const
{
	//DownSample主要从每个fine node计算其对coarse node的贡献
	if (depth == 0) return;
	double cornerValue;
	if (_boundaryType == -1) cornerValue = 0.50;
	else if (_boundaryType == 1) cornerValue = 1.00;
	else                         cornerValue = 0.75;
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
	for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(depth);
#pragma omp parallel for num_threads( threads )
	for (int i = sNodes.nodeCount[depth]; i < sNodes.nodeCount[depth + 1]; i++)
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
		int d, off[3];
		//B样条的coarse函数写成fine函数的线性组合
		UpSampleData usData[3];
		sNodes.treeNodes[i]->depthAndOffset(d, off);
		for (int dd = 0; dd < 3; dd++)
		{
			if (off[dd] == 0) usData[dd] = UpSampleData(1, cornerValue, 0.00);
			else if (off[dd] + 1 == (1 << depth)) usData[dd] = UpSampleData(0, 0.00, cornerValue);
			else if (off[dd] % 2) usData[dd] = UpSampleData(1, 0.75, 0.25);
			else                             usData[dd] = UpSampleData(0, 0.25, 0.75);
		}
		typename TreeOctNode::Neighbors3& neighbors = neighborKey.getNeighbors(sNodes.treeNodes[i]->parent);
		int fine_index = fineConstraints + i - sNodes.nodeCount[depth];
		for (int ii = 0; ii < 2; ii++)
		{
			int _ii = ii + usData[0].start;
			Real cx = Real(usData[0].v[ii]);
			for (int jj = 0; jj < 2; jj++)
			{
				int _jj = jj + usData[1].start;
				Real cxy = Real(cx * usData[1].v[jj]);
				for (int kk = 0; kk < 2; kk++)
				{
					int _kk = kk + usData[2].start;
					TreeOctNode* pNode = neighbors.neighbors[_ii][_jj][_kk];
					if (pNode)
					{
#pragma omp critical
						{
							int coarse_index = coarseConstraints + pNode->nodeData.nodeIndex - sNodes.nodeCount[depth - 1];
							AddSparseVector2(_Q[coarse_index], _Q_index[coarse_index], _Q[fine_index], Real(cxy * usData[2].v[kk]), thre);
						}
						//coarseConstraints[pNode->nodeData.nodeIndex - sNodes.nodeCount[depth - 1]] += C(cxy * usData[2].v[kk]);
					}
				}
			}
		}
	}
}



template< class Real >
void Octree2< Real >::DownSample2List(int depth, const SortedTreeNodes& sNodes, int fineConstraints, int coarseConstraints, std::vector<std::list<MatrixEntry <Real>>>& _Q) const
{
	if (depth == 0) return;
	double cornerValue;
	if (_boundaryType == -1) cornerValue = 0.50;
	else if (_boundaryType == 1) cornerValue = 1.00;
	else                         cornerValue = 0.75;
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
	for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(depth);
#pragma omp parallel for num_threads( threads )
	for (int i = sNodes.nodeCount[depth]; i < sNodes.nodeCount[depth + 1]; i++)
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
		int d, off[3];
		//B样条的coarse函数写成fine函数的线性组合
		UpSampleData usData[3];
		sNodes.treeNodes[i]->depthAndOffset(d, off);
		for (int dd = 0; dd < 3; dd++)
		{
			if (off[dd] == 0) usData[dd] = UpSampleData(1, cornerValue, 0.00);
			else if (off[dd] + 1 == (1 << depth)) usData[dd] = UpSampleData(0, 0.00, cornerValue);
			else if (off[dd] % 2) usData[dd] = UpSampleData(1, 0.75, 0.25);
			else                             usData[dd] = UpSampleData(0, 0.25, 0.75);
		}
		typename TreeOctNode::Neighbors3& neighbors = neighborKey.getNeighbors(sNodes.treeNodes[i]->parent);
		int fine_index = fineConstraints + i - sNodes.nodeCount[depth];
		for (int ii = 0; ii < 2; ii++)
		{
			int _ii = ii + usData[0].start;
			Real cx = Real(usData[0].v[ii]);
			for (int jj = 0; jj < 2; jj++)
			{
				int _jj = jj + usData[1].start;
				Real cxy = Real(cx * usData[1].v[jj]);
				for (int kk = 0; kk < 2; kk++)
				{
					int _kk = kk + usData[2].start;
					TreeOctNode* pNode = neighbors.neighbors[_ii][_jj][_kk];
					if (pNode)
					{
#pragma omp critical
						{
							int coarse_index = coarseConstraints + pNode->nodeData.nodeIndex - sNodes.nodeCount[depth - 1];
							AddSparseList(_Q[coarse_index], _Q[fine_index], Real(cxy * usData[2].v[kk]));
						}
						//coarseConstraints[pNode->nodeData.nodeIndex - sNodes.nodeCount[depth - 1]] += C(cxy * usData[2].v[kk]);
					}
				}
			}
		}
	}
}

template< class Real >
template< class C >
void Octree2< Real >::UpSample(int depth, const SortedTreeNodes& sNodes, ConstPointer(C) coarseCoefficients, Pointer(C) fineCoefficients) const
{
	double cornerValue;
	if (_boundaryType == -1) cornerValue = 0.50;
	else if (_boundaryType == 1) cornerValue = 1.00;
	else                         cornerValue = 0.75;
	if (depth <= _minDepth) return;

	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
	for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(depth - 1);
#pragma omp parallel for num_threads( threads )
	for (int i = sNodes.nodeCount[depth]; i < sNodes.nodeCount[depth + 1]; i++)
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
		bool isInterior = true;
		TreeOctNode* node = sNodes.treeNodes[i];
		int d, off[3];
		UpSampleData usData[3];
		node->depthAndOffset(d, off);
		for (int d = 0; d < 3; d++)
		{
			if (off[d] == 0) usData[d] = UpSampleData(1, cornerValue, 0.00), isInterior = false;
			else if (off[d] + 1 == (1 << depth)) usData[d] = UpSampleData(0, 0.00, cornerValue), isInterior = false;
			else if (off[d] % 2) usData[d] = UpSampleData(1, 0.75, 0.25);
			else                            usData[d] = UpSampleData(0, 0.25, 0.75);
		}
		typename TreeOctNode::Neighbors3& neighbors = neighborKey.getNeighbors(node->parent);
		for (int ii = 0; ii < 2; ii++)
		{
			int _ii = ii + usData[0].start;
			double dx = usData[0].v[ii];
			for (int jj = 0; jj < 2; jj++)
			{
				int _jj = jj + usData[1].start;
				double dxy = dx * usData[1].v[jj];
				for (int kk = 0; kk < 2; kk++)
				{
					int _kk = kk + usData[2].start;
					TreeOctNode* node = neighbors.neighbors[_ii][_jj][_kk];
					if (node)
					{
						double dxyz = dxy * usData[2].v[kk];
						int _i = node->nodeData.nodeIndex;
						fineCoefficients[i - sNodes.nodeCount[depth]] += coarseCoefficients[_i - sNodes.nodeCount[depth - 1]] * Real(dxyz);
					}
				}
			}
		}
	}
}



template< class Real >
void Octree2< Real >::UpSample2(int depth, const SortedTreeNodes& sNodes, int coarseConstraints, int fineConstraints, std::vector<std::vector<MatrixEntry <Real>>>& coefficients) const
{
	double cornerValue;
	if (_boundaryType == -1) cornerValue = 0.50;
	else if (_boundaryType == 1) cornerValue = 1.00;
	else                         cornerValue = 0.75;
	if (depth <= _minDepth) return;

	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
	for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(depth - 1);
#pragma omp parallel for num_threads( threads )
	for (int i = sNodes.nodeCount[depth]; i < sNodes.nodeCount[depth + 1]; i++)
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
		bool isInterior = true;
		TreeOctNode* node = sNodes.treeNodes[i];
		int d, off[3];
		UpSampleData usData[3];
		node->depthAndOffset(d, off);
		for (int d = 0; d < 3; d++)
		{
			if (off[d] == 0) usData[d] = UpSampleData(1, cornerValue, 0.00), isInterior = false;
			else if (off[d] + 1 == (1 << depth)) usData[d] = UpSampleData(0, 0.00, cornerValue), isInterior = false;
			else if (off[d] % 2) usData[d] = UpSampleData(1, 0.75, 0.25);
			else                            usData[d] = UpSampleData(0, 0.25, 0.75);
		}
		typename TreeOctNode::Neighbors3& neighbors = neighborKey.getNeighbors(node->parent);
		for (int ii = 0; ii < 2; ii++)
		{
			int _ii = ii + usData[0].start;
			double dx = usData[0].v[ii];
			for (int jj = 0; jj < 2; jj++)
			{
				int _jj = jj + usData[1].start;
				double dxy = dx * usData[1].v[jj];
				for (int kk = 0; kk < 2; kk++)
				{
					int _kk = kk + usData[2].start;
					TreeOctNode* node = neighbors.neighbors[_ii][_jj][_kk];
					if (node)
					{
						double dxyz = dxy * usData[2].v[kk];
						int _i = node->nodeData.nodeIndex;
						int fine_index = fineConstraints + i - sNodes.nodeCount[depth];
						int coarse_index = coarseConstraints + _i - sNodes.nodeCount[depth - 1];
						AddSparseVector(coefficients[fine_index], coefficients[coarse_index], Real(dxyz));
						//fineCoefficients[i - sNodes.nodeCount[depth]] += coarseCoefficients[_i - sNodes.nodeCount[depth - 1]] * Real(dxyz);
					}
				}
			}
		}
	}
}

template< class Real >
void Octree2< Real >::UpSample3(int depth, const SortedTreeNodes& sNodes, int coarseConstraints, int fineConstraints, std::vector<std::vector<MatrixEntry <Real>>>& coefficients, std::vector < std::unordered_map<int, int>>& coefficients_index, Real thre) const
{
	double cornerValue;
	if (_boundaryType == -1) cornerValue = 0.50;
	else if (_boundaryType == 1) cornerValue = 1.00;
	else                         cornerValue = 0.75;
	if (depth <= _minDepth) return;

	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
	for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(depth - 1);
#pragma omp parallel for num_threads( threads )
	for (int i = sNodes.nodeCount[depth]; i < sNodes.nodeCount[depth + 1]; i++)
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
		bool isInterior = true;
		TreeOctNode* node = sNodes.treeNodes[i];
		int d, off[3];
		UpSampleData usData[3];
		node->depthAndOffset(d, off);
		for (int d = 0; d < 3; d++)
		{
			if (off[d] == 0) usData[d] = UpSampleData(1, cornerValue, 0.00), isInterior = false;
			else if (off[d] + 1 == (1 << depth)) usData[d] = UpSampleData(0, 0.00, cornerValue), isInterior = false;
			else if (off[d] % 2) usData[d] = UpSampleData(1, 0.75, 0.25);
			else                            usData[d] = UpSampleData(0, 0.25, 0.75);
		}
		typename TreeOctNode::Neighbors3& neighbors = neighborKey.getNeighbors(node->parent);
		for (int ii = 0; ii < 2; ii++)
		{
			int _ii = ii + usData[0].start;
			double dx = usData[0].v[ii];
			for (int jj = 0; jj < 2; jj++)
			{
				int _jj = jj + usData[1].start;
				double dxy = dx * usData[1].v[jj];
				for (int kk = 0; kk < 2; kk++)
				{
					int _kk = kk + usData[2].start;
					TreeOctNode* node = neighbors.neighbors[_ii][_jj][_kk];
					if (node)
					{
						double dxyz = dxy * usData[2].v[kk];
						int _i = node->nodeData.nodeIndex;
						int fine_index = fineConstraints + i - sNodes.nodeCount[depth];
						int coarse_index = coarseConstraints + _i - sNodes.nodeCount[depth - 1];
						AddSparseVector2(coefficients[fine_index], coefficients_index[fine_index], coefficients[coarse_index], Real(dxyz), thre);
						//fineCoefficients[i - sNodes.nodeCount[depth]] += coarseCoefficients[_i - sNodes.nodeCount[depth - 1]] * Real(dxyz);
					}
				}
			}
		}
	}
}


template< class Real >
void Octree2< Real >::UpSample2List(int depth, const SortedTreeNodes& sNodes, int coarseConstraints, int fineConstraints, std::vector<std::list<MatrixEntry <Real>>>& coefficients) const
{
	double cornerValue;
	if (_boundaryType == -1) cornerValue = 0.50;
	else if (_boundaryType == 1) cornerValue = 1.00;
	else                         cornerValue = 0.75;
	if (depth <= _minDepth) return;

	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
	for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(depth - 1);
#pragma omp parallel for num_threads( threads )
	for (int i = sNodes.nodeCount[depth]; i < sNodes.nodeCount[depth + 1]; i++)
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
		bool isInterior = true;
		TreeOctNode* node = sNodes.treeNodes[i];
		int d, off[3];
		UpSampleData usData[3];
		node->depthAndOffset(d, off);
		for (int d = 0; d < 3; d++)
		{
			if (off[d] == 0) usData[d] = UpSampleData(1, cornerValue, 0.00), isInterior = false;
			else if (off[d] + 1 == (1 << depth)) usData[d] = UpSampleData(0, 0.00, cornerValue), isInterior = false;
			else if (off[d] % 2) usData[d] = UpSampleData(1, 0.75, 0.25);
			else                            usData[d] = UpSampleData(0, 0.25, 0.75);
		}
		typename TreeOctNode::Neighbors3& neighbors = neighborKey.getNeighbors(node->parent);
		for (int ii = 0; ii < 2; ii++)
		{
			int _ii = ii + usData[0].start;
			double dx = usData[0].v[ii];
			for (int jj = 0; jj < 2; jj++)
			{
				int _jj = jj + usData[1].start;
				double dxy = dx * usData[1].v[jj];
				for (int kk = 0; kk < 2; kk++)
				{
					int _kk = kk + usData[2].start;
					TreeOctNode* node = neighbors.neighbors[_ii][_jj][_kk];
					if (node)
					{
						double dxyz = dxy * usData[2].v[kk];
						int _i = node->nodeData.nodeIndex;
						int fine_index = fineConstraints + i - sNodes.nodeCount[depth];
						int coarse_index = coarseConstraints + _i - sNodes.nodeCount[depth - 1];
						AddSparseList(coefficients[fine_index], coefficients[coarse_index], Real(dxyz), thre);
						//fineCoefficients[i - sNodes.nodeCount[depth]] += coarseCoefficients[_i - sNodes.nodeCount[depth - 1]] * Real(dxyz);
					}
				}
			}
		}
	}
}

// At each point @( depth ), evaluate the met solution @( depth-1 )
template< class Real >
void Octree2< Real >::SetPointValuesFromCoarser(PointInfo& pointInfo, int depth, const SortedTreeNodes& sNodes, ConstPointer(Real) coarseCoefficients)
{
	std::vector< _PointData >& points = pointInfo.points;
	// For every node at the current depth
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
	for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(depth);
#pragma omp parallel for num_threads( threads )
	for (int i = sNodes.nodeCount[depth]; i < sNodes.nodeCount[depth + 1]; i++)
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
		int pIdx = pointInfo.pointIndex(sNodes.treeNodes[i]);
		if (pIdx != -1)
		{
			neighborKey.getNeighbors(sNodes.treeNodes[i]);
			points[pIdx].weightedCoarserValue = _WeightedCoarserFunctionValue(points[pIdx], neighborKey, sNodes.treeNodes[i], coarseCoefficients - _sNodes.nodeCount[depth - 1]);
		}
	}
}
template< class Real >
Real Octree2< Real >::_WeightedCoarserFunctionValue(const _PointData& pointData, const typename TreeOctNode::NeighborKey3& neighborKey, const TreeOctNode* pointNode, ConstPointer(Real) coarseCoefficients) const
{
	double pointValue = 0;
	int depth = pointNode->depth();
	//_boundaryType == -1 && depth == 0
	if (_boundaryType == -1 && depth == 0)
	{
		return Real(-0.5) * pointData.weight;
	}
	if (depth <= _minDepth) return Real(0.);

	Real weight = pointData.weight;
	Point3D< Real > p = pointData.position;

	// Iterate over all basis functions that overlap the point at the coarser resolutions
	{
		int d, _idx[3];
		const typename TreeOctNode::Neighbors3& neighbors = neighborKey.neighbors[depth - 1];
		neighbors.neighbors[1][1][1]->depthAndOffset(d, _idx);
		_idx[0] = BinaryNode::CenterIndex(d, _idx[0] - 1);
		_idx[1] = BinaryNode::CenterIndex(d, _idx[1] - 1);
		_idx[2] = BinaryNode::CenterIndex(d, _idx[2] - 1);

		for (int j = 0; j < 3; j++)
		{
#if ROBERTO_TOLDO_FIX
			double xValue = 0;
			if (_idx[0] + j >= 0 && _idx[0] + j < ((1 << depth) - 1)) xValue = _fData.baseBSplines[_idx[0] + j][2 - j](p[0]);
			else continue;
#else // !ROBERTO_TOLDO_FIX
			double xValue = _fData.baseBSplines[_idx[0] + j][2 - j](p[0]);
#endif // ROBERTO_TOLDO_FIX
			for (int k = 0; k < 3; k++)
			{
#if ROBERTO_TOLDO_FIX
				double xyValue = 0;
				if (_idx[1] + k >= 0 && _idx[1] + k < ((1 << depth) - 1)) xyValue = xValue * _fData.baseBSplines[_idx[1] + k][2 - k](p[1]);
				else continue;
#else // !ROBERTO_TOLDO_FIX
				double xyValue = xValue * _fData.baseBSplines[_idx[1] + k][2 - k](p[1]);
#endif // ROBERTO_TOLDO_FIX
				double _pointValue = 0;
				for (int l = 0; l < 3; l++)
				{
					const TreeOctNode* basisNode = neighbors.neighbors[j][k][l];
#if ROBERTO_TOLDO_FIX
					if (basisNode && basisNode->nodeData.nodeIndex >= 0 && _idx[2] + l >= 0 && _idx[2] + l < ((1 << depth) - 1))
						_pointValue += _fData.baseBSplines[_idx[2] + l][2 - l](p[2]) * double(coarseCoefficients[basisNode->nodeData.nodeIndex]);
#else // !ROBERTO_TOLDO_FIX
					if (basisNode && basisNode->nodeData.nodeIndex >= 0)
						_pointValue += _fData.baseBSplines[_idx[2] + l][2 - l](p[2]) * double(coarseCoefficients[basisNode->nodeData.nodeIndex]);
#endif // ROBERTO_TOLDO_FIX
				}
				pointValue += _pointValue * xyValue;
			}
		}
	}
	//if (_boundaryType == -1)
	if (_boundaryType == -1)
	{
		pointValue -= 0.5;
	}
	return Real(pointValue * weight);
}
template< class Real >
void Octree2< Real >::SetPointConstraintsFromFiner(const PointInfo& pointInfo, int depth, const SortedTreeNodes& sNodes, ConstPointer(Real) finerCoefficients, Pointer(Real) coarserConstraints) const
{
	const std::vector< _PointData >& points = pointInfo.points;
	// Note: We can't iterate over the finer point nodes as the point weights might be
	// scaled incorrectly, due to the adaptive exponent. So instead, we will iterate
	// over the coarser nodes and evaluate the finer solution at the associated points.
	if (!depth) return;
	size_t start = sNodes.nodeCount[depth - 1], end = sNodes.nodeCount[depth], range = end - start;
	memset(coarserConstraints, 0, sizeof(Real) * (sNodes.nodeCount[depth] - sNodes.nodeCount[depth - 1]));
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
	for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(depth - 1);
#pragma omp parallel for num_threads( threads )
	for (int i = sNodes.nodeCount[depth - 1]; i < sNodes.nodeCount[depth]; i++)
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
		int pIdx = pointInfo.pointIndex(sNodes.treeNodes[i]);
		if (pIdx != -1)
		{
			typename TreeOctNode::Neighbors3& neighbors = neighborKey.getNeighbors(sNodes.treeNodes[i]);
			// Evaluate the solution @( depth ) at the current point @( depth-1 )
			{
				Real finerPointValue = _WeightedFinerFunctionValue(points[pIdx], neighborKey, sNodes.treeNodes[i], finerCoefficients - sNodes.nodeCount[depth]);
				Point3D< Real > p = points[pIdx].position;
				// Update constraints for all nodes @( depth-1 ) that overlap the point
				int d, idx[3];
				neighbors.neighbors[1][1][1]->depthAndOffset(d, idx);
				// Set the (offset) index to the top-left-front corner of the 3x3x3 block of b-splines
				// overlapping the point.
				idx[0] = BinaryNode::CenterIndex(d, idx[0] - 1);
				idx[1] = BinaryNode::CenterIndex(d, idx[1] - 1);
				idx[2] = BinaryNode::CenterIndex(d, idx[2] - 1);
				for (int x = 0; x < 3; x++) for (int y = 0; y < 3; y++) for (int z = 0; z < 3; z++)
					if (neighbors.neighbors[x][y][z])
					{
#pragma omp atomic
						coarserConstraints[neighbors.neighbors[x][y][z]->nodeData.nodeIndex - sNodes.nodeCount[depth - 1]] +=
							Real(
								_fData.baseBSplines[idx[0] + x][2 - x](p[0]) *
								_fData.baseBSplines[idx[1] + y][2 - y](p[1]) *
								_fData.baseBSplines[idx[2] + z][2 - z](p[2]) *
								finerPointValue
							);
					}
			}
		}
	}
}
template< class Real >
Real Octree2< Real >::_WeightedFinerFunctionValue(const _PointData& pointData, const typename TreeOctNode::NeighborKey3& neighborKey, const TreeOctNode* pointNode, ConstPointer(Real) finerCoefficients) const
{
	typename TreeOctNode::Neighbors3 childNeighbors;
	double pointValue = 0;
	int depth = pointNode->depth();
	Real weight = pointData.weight;
	Point3D< Real > p = pointData.position;
	neighborKey.getChildNeighbors(p, depth, childNeighbors);
	// Iterate over all finer basis functions that overlap the point at the coarser resolutions
	int d, idx[3];
	{
		Point3D< Real > c;
		Real w;
		neighborKey.neighbors[depth].neighbors[1][1][1]->depthAndOffset(d, idx);
		neighborKey.neighbors[depth].neighbors[1][1][1]->centerAndWidth(c, w);
		d++;
		idx[0] *= 2, idx[1] *= 2, idx[2] *= 2;
		int cIndex = TreeOctNode::CornerIndex(c, p);
		if (cIndex & 1) idx[0]++;
		if (cIndex & 2) idx[1]++;
		if (cIndex & 4) idx[2]++;
	}
	// Center the indexing at the top-left-front corner
	idx[0] = BinaryNode::CenterIndex(d, idx[0] - 1);
	idx[1] = BinaryNode::CenterIndex(d, idx[1] - 1);
	idx[2] = BinaryNode::CenterIndex(d, idx[2] - 1);

	for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++)
	{
		const TreeOctNode* basisNode = childNeighbors.neighbors[j][k][l];
		if (basisNode && basisNode->nodeData.nodeIndex >= 0)
			pointValue +=
			_fData.baseBSplines[idx[0] + j][2 - j](p[0]) *
			_fData.baseBSplines[idx[1] + k][2 - k](p[1]) *
			_fData.baseBSplines[idx[2] + l][2 - l](p[2]) *
			double(finerCoefficients[basisNode->nodeData.nodeIndex]);
	}
	if (_boundaryType == -1) pointValue -= Real(0.5);
	return Real(pointValue * weight);
}
template< class Real >
int Octree2< Real >::GetSliceMatrixAndUpdateConstraints(const PointInfo& pointInfo, SparseMatrix< Real >& matrix, Pointer(Real) constraints, const typename BSplineData< 2 >::Integrator& integrator, int depth, const SortedTreeNodes& sNodes, ConstPointer(Real) metSolution, bool coarseToFine, int nStart, int nEnd)
{
	size_t range = nEnd - nStart;
	Stencil< double, 5 > stencil, stencils[2][2][2];
	SetLaplacianStencil(depth, integrator, stencil);
	SetLaplacianStencils(depth, integrator, stencils);
	matrix.Resize((int)range);
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
	for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(depth);
#pragma omp parallel for num_threads( threads )
	for (int i = 0; i < range; i++)
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
		TreeOctNode* node = sNodes.treeNodes[i + nStart];
		// Get the matrix row size
		bool insetSupported = _boundaryType != 0 || _IsInsetSupported(node);
		typename TreeOctNode::Neighbors5 neighbors5;
		if (insetSupported) neighborKey.getNeighbors(node, neighbors5);
		int count = insetSupported ? GetMatrixRowSize(neighbors5, false) : 1;

		// Allocate memory for the row
#pragma omp critical (matrix_set_row_size)
		{
			matrix.SetRowSize(i, count);
		}

		// Set the row entries
		if (insetSupported) matrix.rowSizes[i] = SetMatrixRowSPR(pointInfo, neighbors5, matrix[i], sNodes.nodeCount[depth], integrator, stencil, false);
		else
		{
			matrix[i][0] = MatrixEntry< Real >(i, Real(1));
			matrix.rowSizes[i] = 1;
		}

		if (depth > _minDepth)
		{
			// Offset the constraints using the solution from lower resolutions.
			int x, y, z, c;
			if (node->parent)
			{
				c = int(node - node->parent->children);
				Cube::FactorCornerIndex(c, x, y, z);
			}
			else x = y = z = 0;
			if (insetSupported && coarseToFine)
			{
				typename TreeOctNode::Neighbors5 pNeighbors5;
				neighborKey.getNeighbors(node->parent, pNeighbors5);
				UpdateConstraintsFromCoarserSPR(pointInfo, neighbors5, pNeighbors5, node, constraints, metSolution, integrator, stencils[x][y][z]);
			}
		}
	}
	return 1;
}
template< class Real >
int Octree2< Real >::GetMatrixAndUpdateConstraints(const PointInfo& pointInfo, SparseSymmetricMatrix< Real >& matrix, SparseSymmetricMatrix< Real >& U, SparseSymmetricMatrix< Real >& H, Vector< Real > &C, const typename BSplineData< 2 >::Integrator& integrator, int depth, const SortedTreeNodes& sNodes, ConstPointer(Real) metSolution, bool coarseToFine)
{
	size_t start = sNodes.nodeCount[depth], end = sNodes.nodeCount[depth + 1], range = end - start;
	Stencil< double, 5 > stencil, stencils[2][2][2];
	SetLaplacianStencil(depth, integrator, stencil);
	SetLaplacianStencils(depth, integrator, stencils);
	matrix.Resize((int)range);
	U.Resize((int)range);
	H.Resize((int)range);
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
	for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(depth);
#pragma omp parallel for num_threads( threads )
	for (int i = 0; i < range; i++)
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
		TreeOctNode* node = sNodes.treeNodes[i + start];
		// Get the matrix row size
		//因为_boundaryType == 0时点云被scale到[0.25, 0.75]，因此_IsInsetSupported用来判断node是否在这个范围内。
		bool insetSupported = _boundaryType != 0 || _IsInsetSupported(node);
		typename TreeOctNode::Neighbors5 neighbors5;
		if (insetSupported) neighborKey.getNeighbors(node, neighbors5);
		int count = insetSupported ? GetMatrixRowSize(neighbors5, true) : 1;

		// Allocate memory for the row
#pragma omp critical (matrix_set_row_size)
		matrix.SetRowSize(i, count);
		U.SetRowSize(i, count);
		H.SetRowSize(i, count);
		// Set the row entries
		int H_size = 0;
		if (insetSupported)
		{
			matrix.rowSizes[i] = SetMatrixRow(pointInfo, neighbors5, matrix[i], U[i], H[i], H_size,  (int)start, integrator, stencil, true);
			U.rowSizes[i] = matrix.rowSizes[i];
		}
		else
		{
			matrix[i][0] = MatrixEntry< Real >(i, Real(1));
			U[i][0] = MatrixEntry< Real >(i, Real(1));
			matrix.rowSizes[i] = 1;
			U.rowSizes[i] = 1;
		}
		if (depth > _minDepth)
		{
			// Offset the constraints using the solution from lower resolutions.
			int x, y, z, c;
			if (node->parent)
			{
				c = int(node - node->parent->children);
				Cube::FactorCornerIndex(c, x, y, z);
			}
			else x = y = z = 0;
			if (insetSupported && coarseToFine)
			{
				typename TreeOctNode::Neighbors5 pNeighbors5;
				neighborKey.getNeighbors(node->parent, pNeighbors5);
				UpdateConstraintsFromCoarser(pointInfo, neighbors5, pNeighbors5, node, C , metSolution, integrator, stencils[x][y][z]);
			}
		}
	}
	return 1;
}

template< class Real >
int Octree2< Real >::GetMatrixAndUpdateConstraints2(const PointInfo& pointInfo, SparseSymmetricMatrix< Real >& matrix, Vector< Real >& C, const typename BSplineData< 2 >::Integrator& integrator, int depth, const SortedTreeNodes& sNodes, ConstPointer(Real) metSolution, bool coarseToFine)
{
	size_t start = sNodes.nodeCount[depth], end = sNodes.nodeCount[depth + 1], range = end - start;
	Stencil< double, 5 > stencil, stencils[2][2][2];
	SetLaplacianStencil(depth, integrator, stencil);
	SetLaplacianStencils(depth, integrator, stencils);
	matrix.Resize((int)range);
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
	for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(depth);
#pragma omp parallel for num_threads( threads )
	for (int i = 0; i < range; i++)
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
		TreeOctNode* node = sNodes.treeNodes[i + start];
		// Get the matrix row size
		//因为_boundaryType == 0时点云被scale到[0.25, 0.75]，因此_IsInsetSupported用来判断node是否在这个范围内。
		bool insetSupported = _boundaryType != 0 || _IsInsetSupported(node);
		typename TreeOctNode::Neighbors5 neighbors5;
		if (insetSupported) neighborKey.getNeighbors(node, neighbors5);
		int count = insetSupported ? GetMatrixRowSize(neighbors5, true) : 1;

		// Allocate memory for the row
#pragma omp critical (matrix_set_row_size)
		matrix.SetRowSize(i, count);
		// Set the row entries
		if (insetSupported)
		{
			matrix.rowSizes[i] = SetMatrixRow2(pointInfo, neighbors5, matrix[i], (int)start, integrator, stencil, true);
		}
		else
		{
			matrix[i][0] = MatrixEntry< Real >(i, Real(1));
			matrix.rowSizes[i] = 1;
		}
		if (depth > _minDepth)
		{
			// Offset the constraints using the solution from lower resolutions.
			int x, y, z, c;
			if (node->parent)
			{
				c = int(node - node->parent->children);
				Cube::FactorCornerIndex(c, x, y, z);
			}
			else x = y = z = 0;
			if (insetSupported && coarseToFine)
			{
				typename TreeOctNode::Neighbors5 pNeighbors5;
				neighborKey.getNeighbors(node->parent, pNeighbors5);
				UpdateConstraintsFromCoarser(pointInfo, neighbors5, pNeighbors5, node, C, metSolution, integrator, stencils[x][y][z]);
			}
		}
	}
	return 1;
}


template< class Real >
int Octree2< Real >::GetMatrixAndUpdateConstraints3_1(const PointInfo& pointInfo, SparseSymmetricMatrix< Real >& matrix, Vector< Real >& C, const typename BSplineData< 2 >::Integrator& integrator, int depth, const SortedTreeNodes& sNodes, ConstPointer(Real) metSolution, bool coarseToFine)
{
	size_t start = sNodes.nodeCount[depth], end = sNodes.nodeCount[depth + 1], range = end - start;
	Stencil< double, 5 > stencil, stencils[2][2][2];
	SetLaplacianStencil(depth, integrator, stencil);
	SetLaplacianStencils(depth, integrator, stencils);
	matrix.Resize((int)range);
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
	for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(depth);
#pragma omp parallel for num_threads( threads )
	for (int i = 0; i < range; i++)
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
		TreeOctNode* node = sNodes.treeNodes[i + start];
		// Get the matrix row size
		//因为_boundaryType == 0时点云被scale到[0.25, 0.75]，因此_IsInsetSupported用来判断node是否在这个范围内。
		bool insetSupported = _boundaryType != 0 || _IsInsetSupported(node);
		typename TreeOctNode::Neighbors5 neighbors5;
		if (insetSupported) neighborKey.getNeighbors(node, neighbors5);
		int count = insetSupported ? GetMatrixRowSize(neighbors5, true) : 1;

		// Allocate memory for the row
#pragma omp critical (matrix_set_row_size)
		matrix.SetRowSize(i, count);
		// Set the row entries
		if (insetSupported)
		{
			matrix.rowSizes[i] = SetMatrixRow2(pointInfo, neighbors5, matrix[i], (int)start, integrator, stencil, true);
		}
		else
		{
			matrix[i][0] = MatrixEntry< Real >(i, Real(1));
			matrix.rowSizes[i] = 1;
		}
	}
	return 1;
}


template< class Real >
int Octree2< Real >::GetMatrixAndUpdateConstraints3_2(const PointInfo& pointInfo, SparseSymmetricMatrix< Real >& matrix, Vector< Real >& C, const typename BSplineData< 2 >::Integrator& integrator, int depth, const SortedTreeNodes& sNodes, ConstPointer(Real) metSolution, bool coarseToFine)
{
	size_t start = sNodes.nodeCount[depth], end = sNodes.nodeCount[depth + 1], range = end - start;
	Stencil< double, 5 > stencil, stencils[2][2][2];
	SetLaplacianStencil(depth, integrator, stencil);
	SetLaplacianStencils(depth, integrator, stencils);
	matrix.Resize((int)range);
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
	for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(depth);
	for (int i = 0; i < range; i++)
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];

		TreeOctNode* node = sNodes.treeNodes[i + start];
		bool insetSupported = _boundaryType != 0 || _IsInsetSupported(node);
		typename TreeOctNode::Neighbors5 neighbors5;
		if (insetSupported) neighborKey.getNeighbors(node, neighbors5);
		if (depth > _minDepth)
		{
			// Offset the constraints using the solution from lower resolutions.
			int x, y, z, c;
			if (node->parent)
			{
				c = int(node - node->parent->children);
				Cube::FactorCornerIndex(c, x, y, z);
			}
			else x = y = z = 0;
			if (insetSupported && coarseToFine)
			{
				typename TreeOctNode::Neighbors5 pNeighbors5;
				neighborKey.getNeighbors(node->parent, pNeighbors5);
				UpdateConstraintsFromFiner(pointInfo, neighbors5, pNeighbors5, node, C, metSolution, integrator, stencils[x][y][z]);
			}
		}
	}
	return 1;
}


template< class Real >
int Octree2< Real >::GetMatrixAndUpdateConstraintsSPR(const PointInfo& pointInfo, SparseSymmetricMatrix< Real >& matrix, Pointer(Real) constraints, const typename BSplineData< 2 >::Integrator& integrator, int depth, const SortedTreeNodes& sNodes, ConstPointer(Real) metSolution, bool coarseToFine)
{
	size_t start = sNodes.nodeCount[depth], end = sNodes.nodeCount[depth + 1], range = end - start;
	Stencil< double, 5 > stencil, stencils[2][2][2];
	//注意B样条是三维的，其grad是(B(x)导,B(y)导,B(z)导)，因此两个grad B的内积是像GetLaplacian那样计算的
	SetLaplacianStencil(depth, integrator, stencil);
	SetLaplacianStencils(depth, integrator, stencils);
	matrix.Resize((int)range);
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
	for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(depth);
#pragma omp parallel for num_threads( threads )
	for (int i = 0; i < range; i++)
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
		TreeOctNode* node = sNodes.treeNodes[i + start];
		// Get the matrix row size
		bool insetSupported = _boundaryType != 0 || _IsInsetSupported(node);
		typename TreeOctNode::Neighbors5 neighbors5;
		if (insetSupported) neighborKey.getNeighbors(node, neighbors5);
		int count = insetSupported ? GetMatrixRowSize(neighbors5, true) : 1;

		// Allocate memory for the row
#pragma omp critical (matrix_set_row_size)
		matrix.SetRowSize(i, count);

		// Set the row entries
		if (insetSupported) matrix.rowSizes[i] = SetMatrixRowSPR(pointInfo, neighbors5, matrix[i], (int)start, integrator, stencil, true);
		else
		{
			matrix[i][0] = MatrixEntry< Real >(i, Real(1));
			matrix.rowSizes[i] = 1;
		}
		if (depth > _minDepth)
		{
			// Offset the constraints using the solution from lower resolutions.
			int x, y, z, c;
			if (node->parent)
			{
				c = int(node - node->parent->children);
				Cube::FactorCornerIndex(c, x, y, z);
			}
			else x = y = z = 0;
			if (insetSupported && coarseToFine)
			{
				typename TreeOctNode::Neighbors5 pNeighbors5;
				neighborKey.getNeighbors(node->parent, pNeighbors5);
				UpdateConstraintsFromCoarserSPR(pointInfo, neighbors5, pNeighbors5, node, constraints, metSolution, integrator, stencils[x][y][z]);
			}
		}
	}
	return 1;
}
template< class Real >
void Octree2< Real >::AddSparseVector(std::vector< MatrixEntry <Real> >& v1, std::vector< MatrixEntry <Real> >& v2) const
{
	for (int i = 0; i < v2.size(); i++)
	{
		bool sameN = false;
		for (int j = 0; j < v1.size(); j++)
		{
			if (v1[j].N == v2[i].N)
			{
				v1[j].Value += v2[i].Value;
				sameN = true;
				break;
			}
		}
		if (sameN == false)
		{
			v1.push_back(v2[i]);
		}
	}
}

template< class Real >
void Octree2< Real >::AddSparseVector2(std::vector< MatrixEntry <Real> >& v1, std::unordered_map<int, int>& v1_index, std::vector< MatrixEntry <Real> >& v2, Real thre) const
{
	for (int i = 0; i < v2.size(); i++)
	{
		if (abs(v2[i].Value) < thre || isnan(v2[i].Value))
		{
			continue;
		}
		std::unordered_map<int, int>::const_iterator got = v1_index.find(v2[i].N);
		if (got == v1_index.end())
		{
			std::pair<int, int> pair(v2[i].N, v1.size());
			v1_index.insert(pair);
			v1.push_back(v2[i]);
		}
		else
		{
			int index = got->second;
			v1[index].Value += v2[i].Value;
		}
	}
}
//如果里面的元素序号都是按顺序的，且没有重复，则有更快速的方法
template< class Real >
void Octree2< Real >::AddSparseList(std::list< MatrixEntry <Real> >& l1, std::list< MatrixEntry <Real> >& l2) const
{
	std::list<MatrixEntry <Real>> ::iterator i1;
	std::list<MatrixEntry <Real>> ::iterator i2;
	i1 = l1.begin();
	i2 = l2.begin();
	//首先,i1,i2都指向两个链表的表头
	while (1)
	{
		//如果l2插入完成，就break
		if (i2 == l2.end())
		{
			break;
		}
		//如果1和2的N相同，把2的value加到1上
		if ((*i1).N == (*i2).N)
		{
			(*i1).Value += (*i2).Value;
			i1++;
			i2++;
		}
		else
		{
			if (i1 == l1.end())
			{
				l1.push_back((*i2));
				i1++;
				i2++;
			}
			else
			{
				if ((*i1).N < (*i2).N)
				{
					i1++;
				}
				else
				{
					l1.insert(i1, *i2);
					i2++;
				}
			}
		}
	}
}


template< class Real >
void Octree2< Real >::AddSparseList(std::list< MatrixEntry <Real> >& l1, std::list< MatrixEntry <Real> >& l2, Real c) const
{
	std::list<MatrixEntry <Real>> ::iterator i1;
	std::list<MatrixEntry <Real>> ::iterator i2;
	i1 = l1.begin();
	i2 = l2.begin();
	//首先,i1,i2都指向两个链表的表头
	while (1)
	{
		//如果l2插入完成，就break
		if (i2 == l2.end())
		{
			break;
		}
		//如果1和2的N相同，把2的value加到1上
		if ((*i1).N == (*i2).N)
		{
			(*i1).Value += (*i2).Value * c;
			i1++;
			i2++;
		}
		else
		{
			if (i1 == l1.end())
			{
				l1.push_back(MatrixEntry <Real>((*i2).N, (*i2).Value * c));
				i1++;
				i2++;
			}
			else
			{
				if ((*i1).N < (*i2).N)
				{
					i1++;
				}
				else
				{
					l1.insert(i1, *i2);
					i2++;
				}
			}
		}
	}
}

template< class Real >
void Octree2< Real >::AddSparseVector(std::vector< MatrixEntry <Real> >& v1, std::vector< MatrixEntry <Real> >& v2, Real c) const
{
	for (int i = 0; i < v2.size(); i++)
	{
		bool sameN = false;
		for (int j = 0; j < v1.size(); j++)
		{
			if (v1[j].N == v2[i].N)
			{
				v1[j].Value += v2[i].Value * c;
				sameN = true;
				break;
			}
		}
		if (sameN == false)
		{
			MatrixEntry <Real> cv2 = v2[i];
			cv2.Value *= c;
			v1.push_back(cv2);
		}
	}
}

template< class Real >
void Octree2< Real >::AddSparseVector2(std::vector< MatrixEntry <Real> >& v1, std::unordered_map<int, int>& v1_index, std::vector< MatrixEntry <Real> >& v2, Real c, Real thre) const
{
	for (int i = 0; i < v2.size(); i++)
	{
		if (abs(v2[i].Value) < thre || isnan(v2[i].Value))
		{
			continue;
		}
		std::unordered_map<int, int>::const_iterator got = v1_index.find(v2[i].N);
		if (got == v1_index.end())
		{
			std::pair<int, int> pair(v2[i].N, v1.size());
			v1_index.insert(pair);
			MatrixEntry <Real> cv2 = v2[i];
			cv2.Value *= c;
			v1.push_back(cv2);
		}
		else
		{
			int index = got->second;
			v1[index].Value += v2[i].Value * c;
		}
	}
}

template< class Real >
std::vector<std::vector<MatrixEntry <Real>>> Octree2< Real >::TransposeSparseMatrix(std::vector<std::vector<MatrixEntry <Real>>> m, int dim2)
{
	std::vector<std::vector<MatrixEntry <Real>>> mT;
	mT.resize(dim2);
	for (int i = 0; i < m.size(); i++)
	{
		for (int j = 0; j < m[i].size(); j++)
		{
			Real value = m[i][j].Value;
			if (value != 0)
			{
				mT[m[i][j].N].push_back(MatrixEntry <Real>(i, value));
			}
		}
	}
	return mT;
}
template< class Real >
Real Octree2< Real >::VectorDot(std::vector<Real> v1, std::vector<Real> v2)
{
	if (v1.size() != v2.size())
	{
		printf("v1 and v2 should have same dimension\n");
		return -1;
	}
	else
	{
		Real dot = Real(0);
#pragma omp parallel for num_threads( threads ) reduction( + : dot )
		for (int i = 0; i < v1.size(); i++)
		{
			if (!isnan(v1[i]) && !isnan(v2[i]))
			{
				dot += (v1[i] * v2[i]);
			}
		}
		return dot;
	}
}

template< class Real >
Real Octree2< Real >::VectorDot(Vector<Real> v1, Vector<Real> v2)
{
	if (v1.Dimensions() != v2.Dimensions())
	{
		printf("v1 and v2 should have same dimension\n");
		return -1;
	}
	else
	{
		Real dot = Real(0);
#pragma omp parallel for num_threads( threads ) reduction( + : dot )
		for (int i = 0; i < v1.Dimensions(); i++)
		{
			dot += (v1[i] * v2[i]);
		}
		return dot;
	}
}


template< class Real >
void Octree2< Real >::LeafMultiplyA(SparseMatrix< Real > Um, SparseMatrix< Real > UmT, SparseMatrix< Real > Qm, SparseMatrix< Real > QmT, SparseMatrix< Real > Pm, SparseMatrix< Real > PmT, Real lambda, Real eta, int threads, Vector<Real> px, Vector<Real> pmu, Vector<Real>& apx, Vector<Real>& apmu)
{
	Vector<Real> apx1, apx2, apx3, apu1, apu2;
	apx1 = UmT.Multiply(Um.Multiply(px, threads), threads);
	apx2 = PmT.Multiply(Pm.Multiply(px, threads), threads);
	apx3 = PmT.Multiply(Qm.Multiply(pmu, threads), threads);

	apu1 = QmT.Multiply(Pm.Multiply(px, threads), threads);
	apu2 = QmT.Multiply(Qm.Multiply(pmu, threads), threads);
#pragma omp parallel for num_threads( threads )
	for (int i = 0; i < apx1.Dimensions(); i++)
	{
		apx1[i] = apx1[i] + eta * apx2[i] - eta * apx3[i];
	}
	for (int i = 0; i < apu1.Dimensions(); i++)
	{
		apu1[i] = -eta * apu1[i] + eta * apu2[i] + lambda * pmu[i];
	}
	apx = apx1;
	apmu = apu1;
}

template< class Real >
void Octree2< Real >::AddOneElement(std::vector<std::vector<MatrixEntry <Real>>> &M, std::vector <std::unordered_map<int, int>> &M_index, int i, int j, Real value)
{
	
	if (i >= M.size())
	{
		return;
	}
	std::unordered_map<int, int>::const_iterator got = M_index[i].find(j);
	if (got == M_index[i].end())
	{
		std::pair<int, int> pair(j, M[i].size());
		M_index[i].insert(pair);
		M[i].push_back(MatrixEntry <Real>(j, value));
	}
	else
	{
		int index = got->second;
		M[i][index].Value += value;
	}
	
	return;
}

template< class Real >
void Octree2< Real >::VectorMatrixToSparseMatrix(std::vector<std::vector<MatrixEntry <Real>>> v, SparseMatrix< Real > &s)
{
	s.Resize(v.size());
	for (int i = 0; i < v.size(); i++)
	{
		s.SetRowSize(i, v[i].size());
		for (int j = 0; j < v[i].size(); j++)
		{
			s[i][j] = v[i][j];
		}
	}
}

//v1 +- alpha * v2
template< class Real >
std::vector<Real> Octree2< Real >::VectorAdd(std::vector<Real> v1, std::vector<Real> v2, bool add, Real alpha)
{
	std::vector<Real> v;
	if (v1.size() == v2.size())
	{
		v.resize(v1.size());
#pragma omp parallel for num_threads( threads )
		for (int i = 0; i < v1.size(); i++)
		{
			if (add)
			{
				v[i] = v1[i] + alpha * v2[i];
			}
			else
			{
				v[i] = v1[i] - alpha * v2[i];
			}
		}
	}
	else
	{
		printf("v1 and v2 should have same dimension\n");
	}
	return v;
}

template< class Real >
Vector<Real> Octree2< Real >::VectorAdd(Vector<Real> v1, Vector<Real> v2, bool add, Real alpha)
{
	Vector<Real> v;
	if (v1.Dimensions() == v2.Dimensions())
	{
		v.Resize(v1.Dimensions());
#pragma omp parallel for num_threads( threads )
		for (int i = 0; i < v1.Dimensions(); i++)
		{
			if (add)
			{
				v[i] = v1[i] + alpha * v2[i];
			}
			else
			{
				v[i] = v1[i] - alpha * v2[i];
			}
		}
	}
	else
	{
		printf("v1 and v2 should have same dimension\n");
	}
	return v;
}

template<class Real>
template<class T2>
Vector<T2> Octree2< Real >::STDToVector(const std::vector<T2> v) const
{
	Vector<T2> V;
	V.Resize(v.size());
	for (int i = 0; i < v.size(); i++)
	{
		V[i] = v[i];
	}
	return V;
}


//Main implementation of the algorithm
template< class Real >
Pointer(Real) Octree2< Real >::SolveSystem2(std::string file_in, PointInfo& pointInfo, std::vector < std::vector< MatrixEntry <Real> >> Qm, std::vector<std::vector<Real>>& coarse_points, std::vector<std::vector<Real>>& coarse_normals, float alpha, float beta, float pointweight, bool oriented_optimized)
{
	//calculate knn regulars
	std::string file_out;
	for (int i = 0; i < file_in.length(); i++)
	{
		if (file_in.length() - i >= 5)
		{
			file_out = file_out + file_in[i];
		}
	}
	file_out = file_out + "_knn.txt";
	int k_calculate = 10;
	int k = 10;
	//Calculating knn by pclknn.exe
	std::string command_knn = ".\\pcl\\pclknn.exe " + file_in + " " + file_out + " " + std::to_string(k) + "\n";
	//std::cout << command_knn << std::endl;
	system(command_knn.c_str());
	std::ifstream ifs;
	ifs.open(file_out);
	if (!ifs) 
	{
		std::cout << "PCL is not installed, or the file path for .\\pcl\\pclknn.exe is wrong." << std::endl;
		exit(1);
	}
	std::vector<std::vector<int>> knn(num_point, std::vector<int>(k, -1));
	for (int i = 0; i < num_point; i++)
	{
		for (int j = 0; j < k; j++)
		{
			ifs >> knn[i][j];
		}
	}
	ifs.close();

	//E_D(n)
	std::vector<std::vector<MatrixEntry <Real>>> Mm(3 * num_point, std::vector<MatrixEntry <Real>>());
	std::vector<std::unordered_map<int, int>>  M_index(3 * num_point);

	//E_COD(n)

	std::vector<std::vector<MatrixEntry <Real>>> Nm(3 * num_point, std::vector<MatrixEntry <Real>>());
	std::vector < std::unordered_map<int, int>>  N_index(3 * num_point);

	Real max_pi_pj_square = 0.00001;
	//rho
	for (int i = 0; i < num_point; i++)
	{
		int index1 = i;
		for (int j = 0; j < k; j++)
		{
			int index2 = knn[i][j];
			std::vector<Real> pi(3, 0);
			pi[0] = all_points[index1][0];
			pi[1] = all_points[index1][1];
			pi[2] = all_points[index1][2];
			std::vector<Real> pj(3, 0);
			pj[0] = all_points[index2][0];
			pj[1] = all_points[index2][1];
			pj[2] = all_points[index2][2];
			Real d = (pi[0] - pj[0]) * (pi[0] - pj[0]) + (pi[1] - pj[1]) * (pi[1] - pj[1]) + (pi[2] - pj[2]) * (pi[2] - pj[2]);
			if (d > max_pi_pj_square)
			{
				max_pi_pj_square = d;
			}
		}
	}
	Real rho = sqrt(max_pi_pj_square) / 2.0;

	//E(D),E(COD)
	for (int i = 0; i < num_point; i++)
	{
		//printf("%d\n", i);
		int index1 = i;
		for (int j = 0; j < k; j++)
		{
			int index2 = knn[i][j];
			std::vector<Real> pi(3, 0);
			pi[0] = all_points[index1][0];
			pi[1] = all_points[index1][1];
			pi[2] = all_points[index1][2];
			std::vector<Real> pj(3, 0);
			pj[0] = all_points[index2][0];
			pj[1] = all_points[index2][1];
			pj[2] = all_points[index2][2];
			Real d = (pi[0] - pj[0]) * (pi[0] - pj[0]) + (pi[1] - pj[1]) * (pi[1] - pj[1]) + (pi[2] - pj[2]) * (pi[2] - pj[2]);
			Real wij = exp(-d / (rho * rho));
			for (int t = 0; t < 3; t++)
			{
				AddOneElement(Mm, M_index, 3 * index1 + t, 3 * index1 + t, wij);
				AddOneElement(Mm, M_index, 3 * index1 + t, 3 * index2 + t, -wij);
				AddOneElement(Mm, M_index, 3 * index2 + t, 3 * index1 + t, -wij);
				AddOneElement(Mm, M_index, 3 * index2 + t, 3 * index2 + t, wij);
			}
			
			std::vector<std::vector<Real>> pipjt(3, std::vector<Real>(3, 0));
			for (int t = 0; t < 3; t++)
			{
				for (int l = 0; l < 3; l++)
				{
					pipjt[t][l] = (pi[t] - pj[t]) * (pi[l] - pj[l]) / d;
					if (!isnan(pipjt[t][l]))
					{
						AddOneElement(Nm, N_index, 3 * index1 + t, 3 * index1 + l, pipjt[t][l] * wij);
						AddOneElement(Nm, N_index, 3 * index1 + t, 3 * index2 + l, pipjt[t][l] * wij);
						AddOneElement(Nm, N_index, 3 * index2 + t, 3 * index1 + l, pipjt[t][l] * wij);
						AddOneElement(Nm, N_index, 3 * index2 + t, 3 * index2 + l, pipjt[t][l] * wij);
					}
				}
			}
		}
	}
	int iter = 0;
	typename BSplineData< 2 >::Integrator integrator;
	_fData.setIntegrator(integrator, _boundaryType == 0);
	Pointer(Real) solution = AllocPointer< Real >(_sNodes.nodeCount[_sNodes.maxDepth]);
	memset(solution, 0, sizeof(Real) * _sNodes.nodeCount[_sNodes.maxDepth]);
	
	std::vector<int> force_zero;
	int force_zero_count = 0;
	
	std::vector<int> zeros_count(_sNodes.nodeCount[_sNodes.maxDepth]);
	std::vector<bool> is_zero(_sNodes.nodeCount[_sNodes.maxDepth]);
	Real boundary_begin = 0.045;
	Real boundary_end = 1.0 - boundary_begin;
	for (int i = 0; i < _sNodes.nodeCount[_sNodes.maxDepth]; i++)
	{
		zeros_count[i] = force_zero_count;
		is_zero[i] = false;
		TreeOctNode* node = _sNodes.treeNodes[i];
		Point3D<Real> center;
		Real width;
		int d, off[3];
		node->depthAndOffset(d, off);
		int res = (1 << d);
		node->centerAndWidth(center, width);
		Real minx = center[0] - 1.5 * width;
		Real maxx = center[0] + 1.5 * width;
		Real miny = center[1] - 1.5 * width;
		Real maxy = center[1] + 1.5 * width;
		Real minz = center[2] - 1.5 * width;
		Real maxz = center[2] + 1.5 * width;
		if (off[0] == 0 || off[0] == (res - 1) || off[1] == 0 || off[1] == (res - 1) || off[2] == 0 || off[2] == (res - 1))
		{
			is_zero[i] = true;
			//printf("%d\n", i);
			force_zero.push_back(i);
			force_zero_count++;
		}
	}
	std::cout << "Establishing the optimization formulation." << std::endl;
	std::vector<int> node_count(_sNodes.maxDepth);
	std::vector<std::vector<MatrixEntry <Real>>> Pm(_sNodes.nodeCount[_sNodes.maxDepth]);
	std::vector < std::unordered_map<int, int>> Pm_index(_sNodes.nodeCount[_sNodes.maxDepth]);
	Real thre = 1e-6;
	for (int d1 = _sNodes.maxDepth - 1; d1 >= _minDepth; d1--)
	{
		if (d1 <= _sNodes.maxDepth - 2)
		{
			DownSample3(d1, _sNodes, _sNodes.nodeCount[d1], _sNodes.nodeCount[d1 - 1], Pm, Pm_index, thre);
		}
		size_t start = _sNodes.nodeCount[d1];
		size_t end = _sNodes.nodeCount[d1 + 1];
		node_count[d1] = end - start;
		Stencil< double, 5 > stencil, stencils[2][2][2];
		SetLaplacianStencil(d1, integrator, stencil);
		SetLaplacianStencils(d1, integrator, stencils);
		std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
		for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(_fData.depth);
		//计算<B^d1, B^{d1 - 1}>
		if (d1 > _minDepth)
		{
#pragma omp parallel for num_threads( threads )
			for (int i = 0; i < node_count[d1]; i++)
			{
				typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
				TreeOctNode* node = _sNodes.treeNodes[i + start];
				bool insetSupported = _boundaryType != 0 || _IsInsetSupported(node);
				int xt, yt, zt, c;
				if (node->parent && insetSupported)
				{
					c = int(node - node->parent->children);
					Cube::FactorCornerIndex(c, xt, yt, zt);
				}
				else
				{
					xt = yt = zt = 0;
					continue;
				}
				typename TreeOctNode::Neighbors5 pNeighbors5;
				neighborKey.getNeighbors(node->parent, pNeighbors5);
				bool isInterior;
				int d, off[3];
				node->depthAndOffset(d, off);
				int o = _boundaryType == 0 ? (1 << (d - 2)) : 0;
				int mn = 4 + o, mx = (1 << d) - 4 - o;
				isInterior = (off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx);
				// Offset the constraints using the solution from lower resolutions.
				int startX = 0, endX = 5, startY = 0, endY = 5, startZ = 0, endZ = 5;
				UpdateCoarserSupportBounds(node, startX, endX, startY, endY, startZ, endZ);

				for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
				{
					if (pNeighbors5.neighbors[x][y][z] && pNeighbors5.neighbors[x][y][z]->nodeData.nodeIndex >= 0)
					{
						const TreeOctNode* _node = pNeighbors5.neighbors[x][y][z];
						std::vector<MatrixEntry< Real >>& P_line = Pm[_node->nodeData.nodeIndex];
						std::unordered_map<int, int>& P_index_line = Pm_index[_node->nodeData.nodeIndex];

						if (isInterior)
						{
#pragma omp critical
							{
								/*
								int _d, _off[3];
								_node->depthAndOffset(_d, _off);
								double ddx = _fData.dot(d, off[0], _d, _off[0], true, true, false);
								double ddy = _fData.dot(d, off[1], _d, _off[1], true, true, false);
								double ddz = _fData.dot(d, off[2], _d, _off[2], true, true, false);
								double vvx = _fData.dot(d, off[0], _d, _off[0], false, false, false);
								double vvy = _fData.dot(d, off[1], _d, _off[1], false, false, false);
								double vvz = _fData.dot(d, off[2], _d, _off[2], false, false, false);

								double dot_value = ddx * vvy * vvz + ddy * vvx * vvz + ddz * vvx * vvy;
								printf("%f, %f\n", dot_value, stencils[xt][yt][zt].values[x][y][z]);
								*/

								P_line.push_back(MatrixEntry< Real >(node->nodeData.nodeIndex, stencils[xt][yt][zt].values[x][y][z]));
								std::pair<int, int> pair(node->nodeData.nodeIndex, P_line.size() - 1);
								P_index_line.insert(pair);
							}
						}
						else
						{
							int _d, _off[3];
							_node->depthAndOffset(_d, _off);
#pragma omp critical
							{
								P_line.push_back(MatrixEntry< Real >(node->nodeData.nodeIndex, GetLaplacian(integrator, d, off, _off, true)));
								std::pair<int, int> pair(node->nodeData.nodeIndex, P_line.size() - 1);
								P_index_line.insert(pair);
							}
						}
						
					}
				}
			}
		}
		//Constructing P (named A in the paper)
		for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(_fData.depth);
#pragma omp parallel for num_threads( threads )
		for (int i = 0; i < node_count[d1]; i++)
		{
			typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
			std::vector<MatrixEntry< Real >>& P_line = Pm[i + start];
			std::unordered_map<int, int>& P_index_line = Pm_index[i + start];
			TreeOctNode* node = _sNodes.treeNodes[i + start];
			// Get the matrix row size
			bool insetSupported = _boundaryType != 0 || _IsInsetSupported(node);
			typename TreeOctNode::Neighbors5 neighbors5;
			if (insetSupported) neighborKey.getNeighbors(node, neighbors5);
			if (insetSupported)
			{
				bool isInterior;
				int d, off[3];
				node->depthAndOffset(d, off);

				int o = _boundaryType == 0 ? (1 << (d - 2)) : 0;
				int mn = 2 + o, mx = (1 << d) - 2 - o;
				isInterior = (off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx);
				int nodeIndex = neighbors5.neighbors[2][2][2]->nodeData.nodeIndex;
				if (!nodeIndex)
				{
					continue;
				}
				Real pointValues[5][5][5];
				if (isInterior) // General case, so try to make fast
				{
					const TreeOctNode* const* _nodes = &neighbors5.neighbors[0][0][0];
					const double* _stencil = &stencil.values[0][0][0];
					Real* _values = &pointValues[0][0][0];
					for (int j = 0; j < 125; j++) _values[j] = Real(_stencil[j]);
#pragma omp critical
					{
						P_line.push_back(MatrixEntry< Real >(nodeIndex, _values[5 * 5 * 2 + 5 * 2 + 2]));
						std::pair<int, int> pair(nodeIndex, P_line.size() - 1);
						P_index_line.insert(pair);
					}
					for (int j = 0; j < 125; j++) if (j != (5 * 5 * 2 + 5 * 2 + 2) && _nodes[j] && _nodes[j]->nodeData.nodeIndex >= nodeIndex)
					{
#pragma omp critical
						{
							P_line.push_back(MatrixEntry< Real >(_nodes[j]->nodeData.nodeIndex, _values[j]));
							std::pair<int, int> pair(_nodes[j]->nodeData.nodeIndex, P_line.size() - 1);
							P_index_line.insert(pair);
						}
					}
				}
				else
				{
					int d, off[3];
					node->depthAndOffset(d, off);
					Real temp = Real(GetLaplacian(integrator, d, off, off, false));
#pragma omp critical
					{
						P_line.push_back(MatrixEntry< Real >(nodeIndex, temp));
						std::pair<int, int> pair(nodeIndex, P_line.size() - 1);
						P_index_line.insert(pair);
					}
					for (int x = 0; x < 5; x++) for (int y = 0; y < 5; y++) for (int z = 0; z < 5; z++)
					{
						if ((x != 2 || y != 2 || z != 2) && neighbors5.neighbors[x][y][z] && neighbors5.neighbors[x][y][z]->nodeData.nodeIndex >= 0 && neighbors5.neighbors[x][y][z]->nodeData.nodeIndex >= nodeIndex)
						{
							const TreeOctNode* _node = neighbors5.neighbors[x][y][z];
							int _d, _off[3];
							_node->depthAndOffset(_d, _off);
							Real temp = Real(GetLaplacian(integrator, d, off, _off, false));
#pragma omp critical
							{
								P_line.push_back(MatrixEntry< Real >(_node->nodeData.nodeIndex, temp));
								std::pair<int, int> pair(_node->nodeData.nodeIndex, P_line.size() - 1);
								P_index_line.insert(pair);
							}
						}
					}
				}
			}
			else
			{
#pragma omp critical
				{
					P_line.push_back(MatrixEntry< Real >(node->nodeData.nodeIndex, 1.0));
					std::pair<int, int> pair(node->nodeData.nodeIndex, P_line.size() - 1);
					P_index_line.insert(pair);
				}
			}
			//P_line.resize(row_size);
		}
	}

	//Constructing U
	std::vector<std::vector<MatrixEntry <Real>>> UmT(_sNodes.nodeCount[_sNodes.maxDepth]);
	std::vector < std::unordered_map<int, int>> UmT_index(_sNodes.nodeCount[_sNodes.maxDepth]);
	std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
	for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(_sNodes.maxDepth - 1);
	std::vector<bool> point_used(num_point, false);
#pragma omp parallel for num_threads( threads )
	for (int i = _sNodes.nodeCount[_sNodes.maxDepth - 1]; i < _sNodes.nodeCount[_sNodes.maxDepth]; i++)
	{
		typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
		TreeOctNode* node = _sNodes.treeNodes[i];
		typename TreeOctNode::Neighbors5 neighbors5;
		//leaf node
		if (node && node->nodeData.leaf_to_pidx.size() != 0)
		{
			bool insetSupported = _boundaryType != 0 || _IsInsetSupported(node);
			if (insetSupported) neighborKey.getNeighbors(node, neighbors5);
		}
		else
		{
			continue;
		}
		/*
		if (node->nodeData.leaf_to_pidx.size() > 1)
		{
			printf("%d\n", node->nodeData.leaf_to_pidx.size());
		}
		*/
		int idx[3]; node->centerIndex(idx);
		Point3D< Real > p;
		for (int pj = 0; pj < node->nodeData.leaf_to_pidx.size(); pj++)
		{
			int p_index = node->nodeData.leaf_to_pidx[pj];
			if (point_used[p_index] == true)
			{
				continue;
			}
			p[0] = all_points[p_index][0];
			p[1] = all_points[p_index][1];
			p[2] = all_points[p_index][2];
			
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						//idx[3] represent the indices of x,y,z coordinates ([0, 2^d - 1])
						const TreeOctNode* _node = neighbors5.neighbors[j + 1][k + 1][l + 1];
						if (_node)
						{
							int _node_index = _node->nodeData.nodeIndex;
							std::vector<MatrixEntry< Real >>& UT_line = UmT[_node_index];
							std::unordered_map<int, int>& UT_index_line = UmT_index[_node_index];
							int _d, _off[3];
							_node->depthAndOffset(_d, _off);
							double valuex = _fData.value(_d, _off[0], 0, p[0], false, false);
							double valuey = _fData.value(_d, _off[1], 0, p[1], false, false);
							double valuez = _fData.value(_d, _off[2], 0, p[2], false, false);
							double value = valuex * valuey * valuez;
#pragma omp critical
							{
								if (abs(value) > thre)
								{
									point_used[p_index] = true;
									UT_line.push_back(MatrixEntry< Real >(p_index, value));
									std::pair<int, int> pair1(p_index, UT_line.size() - 1);
									UT_index_line.insert(pair1);
								}
							}
						}
					}
				}
			}
			TreeOctNode* node_parents = node->parent;
			typename TreeOctNode::Neighbors5 pneighbors5;
			if (node_parents)
			{
				bool insetSupported = _boundaryType != 0 || _IsInsetSupported(node_parents);
				if (insetSupported) neighborKey.getNeighbors(node_parents, pneighbors5);
			}
			else
			{
				continue;
			}
			node_parents->centerIndex(idx);
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						const TreeOctNode* _node = pneighbors5.neighbors[j + 1][k + 1][l + 1];
						if (_node)
						{
							int _node_index = _node->nodeData.nodeIndex;
							std::vector<MatrixEntry< Real >>& UT_line = UmT[_node_index];
							std::unordered_map<int, int>& UT_index_line = UmT_index[_node_index];
							int _d, _off[3];
							_node->depthAndOffset(_d, _off);
							double valuex = _fData.value(_d, _off[0], 0, p[0], false, false);
							double valuey = _fData.value(_d, _off[1], 0, p[1], false, false);
							double valuez = _fData.value(_d, _off[2], 0, p[2], false, false);
							double value = valuex * valuey * valuez;
#pragma omp critical
							{
								if (abs(value) > thre)
								{
									point_used[p_index] = true;
									UT_line.push_back(MatrixEntry< Real >(p_index, value));
									std::pair<int, int> pair1(p_index, UT_line.size() - 1);
									UT_index_line.insert(pair1);
								}
							}
						}
					}
				}
			}
		}
	}
	for (int d1 = _sNodes.maxDepth - 2; d1 >= _minDepth; d1--)
	{
		DownSample3(d1, _sNodes, _sNodes.nodeCount[d1], _sNodes.nodeCount[d1 - 1], UmT, UmT_index, thre);
	}

	std::vector<std::vector<MatrixEntry <Real>>> PmT(_sNodes.nodeCount[_sNodes.maxDepth]);
	PmT = TransposeSparseMatrix(Pm, _sNodes.nodeCount[_sNodes.maxDepth]);
	for (int i = 0; i < _sNodes.nodeCount[_sNodes.maxDepth]; i++)
	{
		for (int j = 0; j < PmT[i].size(); j++)
		{
			if (i != PmT[i][j].N)
			{
				{
					Pm[i].push_back(PmT[i][j]);
				}
			}
		}
	}

	std::vector<int> not_used_treenode;
	for (int i = 0; i < num_point; i++)
	{
		if (point_used[i] == false)
		{
			not_used_treenode.push_back(i);
		}
	}
#pragma omp parallel for num_threads( threads )
	for (int i = 0; i < _sNodes.nodeCount[_sNodes.maxDepth]; i++)
	{
		TreeOctNode* node = _sNodes.treeNodes[i];
		int d, off[3];
		node->depthAndOffset(d, off);
		int node_index = node->nodeData.nodeIndex;

		std::vector<MatrixEntry< Real >>& UT_line = UmT[node_index];
		std::unordered_map<int, int>& UT_index_line = UmT_index[node_index];
		Point3D<Real> center;
		Real width;
		node->centerAndWidth(center, width);
		Real minx = center[0] - 1.5 * width;
		Real maxx = center[0] + 1.5 * width;
		Real miny = center[1] - 1.5 * width;
		Real maxy = center[1] + 1.5 * width;
		Real minz = center[2] - 1.5 * width;
		Real maxz = center[2] + 1.5 * width;
		for (int j = 0; j < not_used_treenode.size(); j++)
		{
			Point3D<Real> p = all_points[not_used_treenode[j]];
			if (p[0] >= maxx || p[0] <= minx || p[1] >= maxy || p[1] <= miny || p[2] >= maxz || p[2] <= minz)
			{
				continue;
			}
			double valuex = _fData.value(d, off[0], 0, p[0], false, false);
			double valuey = _fData.value(d, off[1], 0, p[1], false, false);
			double valuez = _fData.value(d, off[2], 0, p[2], false, false);
			double value = valuex * valuey * valuez;
#pragma omp critical
			{
				if (abs(value) > thre)
				{
					UT_line.push_back(MatrixEntry< Real >(not_used_treenode[j], value));
					std::pair<int, int> pair1(not_used_treenode[j], UT_line.size() - 1);
					UT_index_line.insert(pair1);
					//printf("value: %d, %f\n", not_used_treenode[j], value);
				}
			}
		}
	}
	
	SparseMatrix< Real > UT;
	SparseMatrix< Real > Q;
	SparseMatrix< Real > P;
	//Data structure changes from vector to SparseMatrix.
	VectorMatrixToSparseMatrix(Pm, P);
	VectorMatrixToSparseMatrix(UmT, UT);
	VectorMatrixToSparseMatrix(Qm, Q);
	std::vector<Real> d(num_point, Real(0.5));

	std::vector< std::vector<MatrixEntry <Real> >> Um, QmT;
	SparseMatrix< Real > U, QT;

	Um = TransposeSparseMatrix(UmT, num_point);
	QmT = TransposeSparseMatrix(Qm, 3 * num_point);
	VectorMatrixToSparseMatrix(Um, U);
	VectorMatrixToSparseMatrix(QmT, QT);

	SparseMatrix< Real > M, N;
	VectorMatrixToSparseMatrix(Mm, M);
	VectorMatrixToSparseMatrix(Nm, N);
	int x_dim = _sNodes.nodeCount[_sNodes.maxDepth];
	int u_dim = 3 * num_point;
	int iter_cg = 300;
	//parameters
	Real eta = alpha;
	Real lambda = beta;
	Real lambda2 = beta / 2.0;
	Real lambda3 = beta;
	
	int threads = 10;
	Vector<Real> dv = STDToVector(d);
	Vector<Real> bx = UT.Multiply(dv, threads);
	std::vector<Real> bmu0(3 * num_point, Real(0.0));
	Vector<Real> bu = STDToVector(bmu0);
	Vector<Real> x = STDToVector(std::vector<Real>(x_dim, Real(0.001)));
	Vector<Real> u = STDToVector(std::vector<Real>(u_dim, Real(0.0)));
	//CG
	Vector<Real> rx, ru, px, pu, ax, au;
	Vector<Real> apx1, apx2, apx3, apu1, apu2, apu3, apu4;
	px = x;
	pu = u;
	for (int j = 0; j < force_zero_count; j++)
	{
		bx[force_zero[j]] = 0;
		px[force_zero[j]] = 0;
	}
	Vector<Real> Ux, Px, Qu, Mu, Nu;
	Ux = U.Multiply(px, threads);
	Px = P.Multiply(px, threads);
	Qu = Q.Multiply(pu, threads);
	Mu = M.Multiply(pu, threads);
	Nu = N.Multiply(pu, threads);
	apx1 = UT.Multiply(Ux, threads);
	apx2 = P.Multiply(Px, threads);
	apx3 = P.Multiply(Qu, threads);
	//printf("dim: %d, %d, %d, %d\n", Px.Dimensions(), apx1.Dimensions(), apx2.Dimensions(), apx3.Dimensions());
	apu1 = QT.Multiply(Px, threads);
	apu2 = QT.Multiply(Qu, threads);
	apu3 = Mu;
	apu4 = Nu;
	ax = apx1 + apx2 * eta - apx3 * eta;
	au = apu1 * (-eta) + apu2 * eta + pu * lambda + apu3 * lambda2 + apu4 * lambda3;
	for (int j = 0; j < force_zero_count; j++)
	{
		ax[force_zero[j]] = 0;
	}
	rx = bx - ax;
	ru = bu - au;
	px = rx;
	pu = ru;
	Real alphak, betak, alphakm, betakm, alphakd, betakd;
	std::cout << "CG optimization." << std::endl;
	for (int i = 0; i < iter_cg; i++)
	{
		
		alphakm = VectorDot(rx, rx) + VectorDot(ru, ru);
		for (int j = 0; j < force_zero_count; j++)
		{
			px[force_zero[j]] = 0;
		}
		Ux = U.Multiply(px, threads);
		Px = P.Multiply(px, threads);
		Qu = Q.Multiply(pu, threads);
		Mu = M.Multiply(pu, threads);
		Nu = N.Multiply(pu, threads);
		
		alphakd = VectorDot(Ux, Ux) + eta * VectorDot(Px, Px) + eta * VectorDot(Qu, Qu)
			+ lambda * VectorDot(pu, pu) - 2 * eta * VectorDot(Px, Qu)
			+ lambda2 * VectorDot(pu, Mu) + lambda3 * VectorDot(pu, Nu);
		if (isnan(alphakd) || isnan(alphakm))
		{
			std::cout << "Detecting nan values, please try again." << std::endl;
			exit(1);
		}
		if (alphakd < Real(0.0001))
		{
			alphakd = Real(0.0001);
		}
		alphak = alphakm / alphakd;
		

		x += (px * alphak);
		u += (pu * alphak);
		Vector<Real> apx1, apx2, apx3, apu1, apu2;
		apx1 = UT.Multiply(Ux, threads);
		apx2 = P.Multiply(Px, threads);
		apx3 = P.Multiply(Qu, threads);

		apu1 = QT.Multiply(Px, threads);
		apu2 = QT.Multiply(Qu, threads);
		apu3 = Mu;
		apu4 = Nu;
		ax = apx1 + apx2 * eta - apx3 * eta;
		au = apu1 * (-eta) + apu2 * eta + pu * lambda + apu3 * lambda2 + apu4 * lambda3;
		for (int j = 0; j < force_zero_count; j++)
		{
			ax[force_zero[j]] = 0;
		}
		rx -= (ax * alphak);
		ru -= (au * alphak);

		betakm = VectorDot(rx, rx) + VectorDot(ru, ru);

		//printf("residual: %f\n", betakm);
		betakd = alphakm;
		if (betakd < Real(0.0001))
		{
			betakd = Real(0.0001);
		}
		betak = betakm / betakd;
		px = rx + px * betak;
		pu = ru + pu * betak;

		/*
		Vector<Real> px_qu = P.Multiply(x, threads) - Q.Multiply(u, threads);
		Real dotpxqu = VectorDot(px_qu, px_qu);
		printf("pxqu:%f\n", eta * dotpxqu);

		Vector<Real> ux_d = U.Multiply(x, threads) - dv;
		Real dotuxd = VectorDot(ux_d, ux_d);
		printf("dotuxd:%f\n", dotuxd);

		Real dotuu = VectorDot(u, u);
		printf("dotuu:%f\n", lambda * dotuu);

		Real dotMu = VectorDot(u, N.Multiply(u, threads));
		printf("dotMu:%f\n", lambda2 * dotMu);

		Real dotNu = VectorDot(u, M.Multiply(u, threads));
		printf("dotNu:%f\n", lambda3 * dotNu);
		*/
	}
	coarse_points.resize(num_point);
	coarse_normals.resize(num_point);
	for (int i = 0; i < x.Dimensions(); i++)
	{
		solution[i] = x[i];
	}
	for (int j = 0; j < num_point; j++)
	{
		coarse_points[j].resize(3);
		coarse_normals[j].resize(3);
		Real l = sqrt(u[3 * j] * u[3 * j] + u[3 * j + 1] * u[3 * j + 1] + u[3 * j + 2] * u[3 * j + 2]);
		//printf("%f\n", l);
		if (l > 1e-6)
		{
			u[3 * j] = u[3 * j] / l;
			u[3 * j + 1] = u[3 * j + 1] / l;
			u[3 * j + 2] = u[3 * j + 2] / l;
		}
		Real ori_x = all_points[j][0] * _scale + _center[0];
		Real ori_y = all_points[j][1] * _scale + _center[1];
		Real ori_z = all_points[j][2] * _scale + _center[2];
		coarse_points[j][0] = ori_x;
		coarse_points[j][1] = ori_y;
		coarse_points[j][2] = ori_z;
		Real dot = u[3 * j] * all_normals[j][0] + u[3 * j + 1] * all_normals[j][1] + u[3 * j + 2] * all_normals[j][2];
		if (oriented_optimized == false)
		{
			//*-1 for showing outward normals.
			coarse_normals[j][0] = -u[3 * j];
			coarse_normals[j][1] = -u[3 * j + 1];
			coarse_normals[j][2] = -u[3 * j + 2];

		}
		else
		{
			//*-1 for showing outward normals.
			if (dot >= 0)
			{
				coarse_normals[j][0] = -all_normals[j][0];
				coarse_normals[j][1] = -all_normals[j][1];
				coarse_normals[j][2] = -all_normals[j][2];
			}
			else
			{
				coarse_normals[j][0] = all_normals[j][0];
				coarse_normals[j][1] = all_normals[j][1];
				coarse_normals[j][2] = all_normals[j][2];
			}
		}
	}
	return solution;
}

template< class Real >
void Octree2< Real >::SolveSPR(PointInfo& pointInfo, Vector<Real> constrains, Vector<Real>& solutions, bool showResidual, int iters, int maxSolveDepth, int cgDepth = 0, double cgAccuracy = 0)
{
	
	int iter = 0;
	typename BSplineData< 2 >::Integrator integrator;
	_fData.setIntegrator(integrator, _boundaryType == 0);
	iters = std::max< int >(0, iters);
	if (_boundaryType == 0) maxSolveDepth++, cgDepth++;

	Pointer(Real) solution = AllocPointer< Real >(_sNodes.nodeCount[_sNodes.maxDepth]);
	memset(solution, 0, sizeof(Real) * _sNodes.nodeCount[_sNodes.maxDepth]);

	solution[0] = 0;

	std::vector< Real > metSolution(_sNodes.nodeCount[_sNodes.maxDepth - 1], 0);

	Pointer(Real) constraints = AllocPointer< Real >(_sNodes.nodeCount[_sNodes.maxDepth]);
	for (int i = 0; i < _sNodes.nodeCount[_sNodes.maxDepth]; i++)
	{
		constraints[i] = constrains[i];
	}
	for (int d = _minDepth; d < _sNodes.maxDepth; d++)
	{
		//DumpOutput("Depth[%d/%d]: %d\n", _boundaryType == 0 ? d - 1 : d, _boundaryType == 0 ? _sNodes.maxDepth - 2 : _sNodes.maxDepth - 1, _sNodes.nodeCount[d + 1] - _sNodes.nodeCount[d]);
		if (d == _minDepth)
			_SolveSPRCG(pointInfo, d, integrator, _sNodes, solution, constraints, GetPointer(metSolution), _sNodes.nodeCount[_minDepth + 1] - _sNodes.nodeCount[_minDepth], true, showResidual, NULL, NULL, NULL);
		else
		{
			if (d > cgDepth) iter += _SolveSPRGS(pointInfo, d, integrator, _sNodes, solution, constraints, GetPointer(metSolution), d > maxSolveDepth ? 0 : iters, true, showResidual, NULL, NULL, NULL);
			else            iter += _SolveSPRCG(pointInfo, d, integrator, _sNodes, solution, constraints, GetPointer(metSolution), d > maxSolveDepth ? 0 : iters, true, showResidual, NULL, NULL, NULL, cgAccuracy);
		}
	}
	for (int i = 0; i < _sNodes.nodeCount[_sNodes.maxDepth]; i++)
	{
		solutions[i] = solution[i];
	}
}

template< class Real >
void Octree2< Real >::_setMultiColorIndices(int start, int end, std::vector< std::vector< int > >& indices) const
{
	const int modulus = 3;
	indices.resize(modulus * modulus * modulus);
	int count[modulus * modulus * modulus];
	memset(count, 0, sizeof(int) * modulus * modulus * modulus);
#pragma omp parallel for num_threads( threads )
	for (int i = start; i < end; i++)
	{
		int d, off[3];
		_sNodes.treeNodes[i]->depthAndOffset(d, off);
		int idx = (modulus * modulus) * (off[2] % modulus) + modulus * (off[1] % modulus) + (off[0] % modulus);
#pragma omp atomic
		count[idx]++;
	}

	for (int i = 0; i < modulus * modulus * modulus; i++) indices[i].reserve(count[i]), count[i] = 0;

	for (int i = start; i < end; i++)
	{
		int d, off[3];
		_sNodes.treeNodes[i]->depthAndOffset(d, off);
		int idx = (modulus * modulus) * (off[2] % modulus) + modulus * (off[1] % modulus) + (off[0] % modulus);
		indices[idx].push_back(_sNodes.treeNodes[i]->nodeData.nodeIndex - start);
	}
}
template< class Real >
int Octree2< Real >::_SolveSystemGS(PointInfo& pointInfo, int depth, const typename BSplineData< 2 >::Integrator& integrator, const SortedTreeNodes& sNodes, Pointer(Real) solution, Pointer(Real) constraints, Pointer(Real) metSolutionConstraints, int iters, bool coarseToFine, bool showResidual, double* bNorm2, double* inRNorm2, double* outRNorm2, bool forceSilent)
{
	Pointer(Real) metSolution = NullPointer< Real >();
	Pointer(Real) metConstraints = NullPointer< Real >();
	if (coarseToFine) metSolution = metSolutionConstraints;	// This stores the up-sampled solution up to depth-2
	else               metConstraints = metSolutionConstraints; // This stores the down-sampled constraints up to depth

	double _maxMemoryUsage = maxMemoryUsage;
	maxMemoryUsage = 0;
	Vector< Real > X, B;
	int slices = 1 << depth;
	double systemTime = 0., solveTime = 0., updateTime = 0., evaluateTime = 0.;
	std::vector< int > offsets(slices + 1, 0);
	for (int i = sNodes.nodeCount[depth]; i < sNodes.nodeCount[depth + 1]; i++)
	{
		int d, off[3];
		sNodes.treeNodes[i]->depthAndOffset(d, off);
		offsets[off[2]]++;
	}
	for (int i = 1; i < slices; i++)  offsets[i] += offsets[i - 1];
	for (int i = slices; i >= 1; i--) offsets[i] = offsets[i - 1];
	offsets[0] = 0;

	X.Resize(sNodes.nodeCount[depth + 1] - sNodes.nodeCount[depth]);
	B.Resize(sNodes.nodeCount[depth + 1] - sNodes.nodeCount[depth]);
	if (coarseToFine)
	{
		if (depth > _minDepth)
		{
			// Up-sample the cumulative change in solution @(depth-2) into the cumulative change in solution @(depth-1)
			if (depth - 2 >= _minDepth) UpSample(depth - 1, sNodes, (ConstPointer(Real))metSolution + _sNodes.nodeCount[depth - 2], metSolution + _sNodes.nodeCount[depth - 1]);
			// Add in the change in solution @(depth-1)
#pragma omp parallel for num_threads( threads )
			for (int i = _sNodes.nodeCount[depth - 1]; i < _sNodes.nodeCount[depth]; i++) metSolution[i] += solution[i];
			// Evaluate the points @(depth) using the cumulative change in solution @(depth-1)
			if (_constrainValues)
			{
				evaluateTime = Time();
				SetPointValuesFromCoarser(pointInfo, depth, sNodes, metSolution + _sNodes.nodeCount[depth - 1]);
				evaluateTime = Time() - evaluateTime;
			}
		}
	}
	else if (depth < _sNodes.maxDepth - 1)
		for (int i = _sNodes.nodeCount[depth]; i < _sNodes.nodeCount[depth + 1]; i++) constraints[i] -= metConstraints[i];
	// Initialize with the previously computed solution
#pragma omp parallel for num_threads( threads )
	for (int i = _sNodes.nodeCount[depth]; i < _sNodes.nodeCount[depth + 1]; i++) X[i - _sNodes.nodeCount[depth]] = solution[i];
	double bNorm = 0, inRNorm = 0, outRNorm = 0;
	if (depth >= _minDepth)
	{
		int frontOffset = (showResidual || inRNorm2) ? 2 : 0;
		int backOffset = (showResidual || outRNorm2) ? 2 : 0;
		int solveSlices = std::min< int >(2 * iters - 1, slices), matrixSlices = std::max< int >(1, std::min< int >(solveSlices + frontOffset + backOffset, slices));
		std::vector< SparseMatrix< Real > > _M(matrixSlices);
		std::vector< std::vector< std::vector< int > > > __mcIndices(std::max< int >(0, solveSlices));

		int dir = coarseToFine ? -1 : 1, start = coarseToFine ? slices - 1 : 0, end = coarseToFine ? -1 : slices;
		for (int frontSlice = start - frontOffset * dir, backSlice = frontSlice - 2 * (iters - 1) * dir; backSlice != end + backOffset * dir; frontSlice += dir, backSlice += dir)
		{
			double t;
			if (frontSlice + frontOffset * dir >= 0 && frontSlice + frontOffset * dir < slices)
			{
				int s = frontSlice + frontOffset * dir, _s = s % matrixSlices;
				t = Time();
				GetSliceMatrixAndUpdateConstraints(pointInfo, _M[_s], constraints, integrator, depth, sNodes, metSolution, coarseToFine, offsets[s] + sNodes.nodeCount[depth], offsets[s + 1] + sNodes.nodeCount[depth]);
				systemTime += Time() - t;
				Pointer(TreeOctNode*) const nodes = sNodes.treeNodes + sNodes.nodeCount[depth];
				for (int i = offsets[s]; i < offsets[s + 1]; i++)
				{
					if (_boundaryType != 0 || _IsInsetSupported(nodes[i])) B[i] = constraints[nodes[i]->nodeData.nodeIndex];
					else                                                    B[i] = Real(0);
				}
				if (showResidual || inRNorm2)
#pragma omp parallel for num_threads( threads ) reduction( + : bNorm , inRNorm )
					for (int j = 0; j < _M[_s].rows; j++)
					{
						Real temp = Real(0);
						ConstPointer(MatrixEntry< Real >) start = _M[_s][j];
						ConstPointer(MatrixEntry< Real >) end = start + _M[_s].rowSizes[j];
						ConstPointer(MatrixEntry< Real >) e;
						for (e = start; e != end; e++) temp += X[e->N] * e->Value;
						Real b = B[j + offsets[s]];
						bNorm += b * b;
						inRNorm += (temp - b) * (temp - b);
					}
				else if (bNorm2)
#pragma omp parallel for num_threads( threads ) reduction( + : bNorm )
					for (int j = 0; j < _M[_s].rows; j++)
					{
						Real b = B[j + offsets[s]];
						bNorm += b * b;
					}
			}
			t = Time();
			if (iters && frontSlice >= 0 && frontSlice < slices)
			{
				int s = frontSlice, _s = s % matrixSlices, __s = s % solveSlices;
				for (int i = 0; i<int(__mcIndices[__s].size()); i++) __mcIndices[__s][i].clear();
				_setMultiColorIndices(sNodes.nodeCount[depth] + offsets[s], sNodes.nodeCount[depth] + offsets[s + 1], __mcIndices[__s]);
			}
			for (int slice = frontSlice; slice * dir >= backSlice * dir; slice -= 2 * dir)
				if (slice >= 0 && slice < slices)
				{
					int s = slice, _s = s % matrixSlices, __s = s % solveSlices;
					SparseMatrix< Real >::SolveGS(__mcIndices[__s], _M[_s], B, X, !coarseToFine, threads, offsets[s]);
				}
			solveTime += Time() - t;
			if ((showResidual || outRNorm2) && backSlice - backOffset * dir >= 0 && backSlice - backOffset * dir < slices)
			{
				int s = backSlice - backOffset * dir, _s = s % matrixSlices;
#pragma omp parallel for num_threads( threads ) reduction( + : outRNorm )
				for (int j = 0; j < _M[_s].rows; j++)
				{
					Real temp = Real(0);
					ConstPointer(MatrixEntry< Real >) start = _M[_s][j];
					ConstPointer(MatrixEntry< Real >) end = start + _M[_s].rowSizes[j];
					ConstPointer(MatrixEntry< Real >) e;
					for (e = start; e != end; e++) temp += X[e->N] * e->Value;
					Real b = B[j + offsets[s]];
					outRNorm += (temp - b) * (temp - b);
				}
			}
		}
	}

	if (bNorm2) bNorm2[depth] = bNorm;
	if (inRNorm2) inRNorm2[depth] = inRNorm;
	if (outRNorm2) outRNorm2[depth] = outRNorm;
	if (showResidual && iters)
	{
		for (int i = 0; i < depth; i++) printf("  ");
		printf("GS: %.4e -> %.4e -> %.4e (%.2e) [%d]\n", sqrt(bNorm), sqrt(inRNorm), sqrt(outRNorm), sqrt(outRNorm / bNorm), iters);
	}

	// Copy the old solution into the buffer, write in the new solution, compute the change, and update the met constraints
#pragma omp parallel for num_threads( threads )
	for (int i = sNodes.nodeCount[depth]; i < sNodes.nodeCount[depth + 1]; i++) solution[i] = X[i - sNodes.nodeCount[depth]];
	if (!coarseToFine && depth > _minDepth)
	{
		// Explicitly compute the restriction of the met solution onto the coarser nodes
		// and down-sample the previous accumulation
		{
			UpdateConstraintsFromFiner(integrator, depth, sNodes, GetPointer(X), metConstraints + sNodes.nodeCount[depth - 1]);
			if (_constrainValues) SetPointConstraintsFromFiner(pointInfo, depth, sNodes, GetPointer(X), metConstraints + sNodes.nodeCount[depth - 1]);
			if (depth < sNodes.maxDepth - 1) DownSample(depth, sNodes, (ConstPointer(Real))metConstraints + sNodes.nodeCount[depth], metConstraints + sNodes.nodeCount[depth - 1]);
		}
	}

	MemoryUsage();
	if (!forceSilent) DumpOutput("\tEvaluated / Got / Solved in: %6.3f / %6.3f / %6.3f\t(%.3f MB)\n", evaluateTime, systemTime, solveTime, float(maxMemoryUsage));
	maxMemoryUsage = std::max< double >(maxMemoryUsage, _maxMemoryUsage);

	return iters;
}


template< class Real >
int Octree2< Real >::_SolveSPRCG(PointInfo& pointInfo, int depth, const typename BSplineData< 2 >::Integrator& integrator, const SortedTreeNodes& sNodes, Pointer(Real) solution, Pointer(Real) constraints, Pointer(Real) metSolutionConstraints, int iters, bool coarseToFine, bool showResidual, double* bNorm2, double* inRNorm2, double* outRNorm2, double accuracy)
{
	Pointer(Real) metSolution = NullPointer< Real >();
	Pointer(Real) metConstraints = NullPointer< Real >();
	if (coarseToFine) metSolution = metSolutionConstraints;	// This stores the up-sampled solution up to depth-2
	else               metConstraints = metSolutionConstraints; // This stores the down-sampled constraints up to depth
	double _maxMemoryUsage = maxMemoryUsage;
	maxMemoryUsage = 0;
	int iter = 0;
	Vector< Real > X, B;
	SparseSymmetricMatrix< Real > M;
	double systemTime = 0., solveTime = 0., updateTime = 0., evaluateTime = 0.;
	X.Resize(sNodes.nodeCount[depth + 1] - sNodes.nodeCount[depth]);
	if (coarseToFine)
	{
		if (depth > _minDepth)
		{
			// Up-sample the cumulative change in solution @(depth-2) into the cumulative change in solution @(depth-1)
			if (depth - 2 >= _minDepth) UpSample(depth - 1, sNodes, (ConstPointer(Real))metSolution + _sNodes.nodeCount[depth - 2], metSolution + _sNodes.nodeCount[depth - 1]);
			// Add in the change in solution @(depth-1)
#pragma omp parallel for num_threads( threads )
			for (int i = _sNodes.nodeCount[depth - 1]; i < _sNodes.nodeCount[depth]; i++) metSolution[i] += solution[i];
			// Evaluate the points @(depth) using the cumulative change in solution @(depth-1)
			/*
			if (_constrainValues)
			{
				evaluateTime = Time();
				SetPointValuesFromCoarser(pointInfo, depth, sNodes, metSolution + _sNodes.nodeCount[depth - 1]);
				evaluateTime = Time() - evaluateTime;
			}
			*/
		}
	}
	else if (depth < _sNodes.maxDepth - 1)
		for (int i = _sNodes.nodeCount[depth]; i < _sNodes.nodeCount[depth + 1]; i++) constraints[i] -= metConstraints[i];  //coarse to fine的时候这就是metSolution \bat{x}^{d-1}
	// Initialize with the previously computed solution
#pragma omp parallel for num_threads( threads )
	for (int i = _sNodes.nodeCount[depth]; i < _sNodes.nodeCount[depth + 1]; i++) X[i - _sNodes.nodeCount[depth]] = solution[i];
	systemTime = Time();
	{
		// Get the system matrix (and adjust the right-hand-side based on the coarser solution if prolonging)
		if (coarseToFine) GetMatrixAndUpdateConstraintsSPR(pointInfo, M, constraints, integrator, depth, sNodes, metSolution, true);
		else               GetMatrixAndUpdateConstraintsSPR(pointInfo, M, constraints, integrator, depth, sNodes, NullPointer< Real >(), false);
		// Set the constraint vector
		B.Resize(sNodes.nodeCount[depth + 1] - sNodes.nodeCount[depth]);
		for (int i = sNodes.nodeCount[depth]; i < sNodes.nodeCount[depth + 1]; i++)
			if (_boundaryType != 0 || _IsInsetSupported(sNodes.treeNodes[i])) B[i - sNodes.nodeCount[depth]] = constraints[i];
			else                                                               B[i - sNodes.nodeCount[depth]] = Real(0);
	}
	systemTime = Time() - systemTime;

	solveTime = Time();
	// Solve the linear system
	accuracy = Real(accuracy / 100000) * M.rows;
	int res = 1 << depth;

	MapReduceVector< Real > mrVector;
	mrVector.resize(threads, M.rows);
	bool addDCTerm = (M.rows == res * res * res && !_constrainValues && _boundaryType != -1);
	double bNorm, inRNorm, outRNorm;
	if (showResidual || bNorm2) bNorm = B.Norm(2);
	if (showResidual || inRNorm2) inRNorm = (addDCTerm ? (B - M * X - X.Average()) : (B - M * X)).Norm(2);

	if (_boundaryType == 0 && depth > 3) res -= 1 << (depth - 2);
	if (iters) iter += SparseSymmetricMatrix< Real >::SolveCG(M, B, iters, X, mrVector, Real(accuracy), 0, addDCTerm);
	solveTime = Time() - solveTime;
	if (showResidual || outRNorm2) outRNorm = (addDCTerm ? (B - M * X - X.Average()) : (B - M * X)).Norm(2);
	if (bNorm2) bNorm2[depth] = bNorm * bNorm;
	if (inRNorm2) inRNorm2[depth] = inRNorm * inRNorm;
	if (outRNorm2) outRNorm2[depth] = outRNorm * outRNorm;
	if (showResidual && iters)
	{
		for (int i = 0; i < depth; i++) printf("  ");
		printf("CG: %.4e -> %.4e -> %.4e (%.2e) [%d]\n", bNorm, inRNorm, outRNorm, outRNorm / bNorm, iter);
	}

	// Copy the old solution into the buffer, write in the new solution, compute the change, and update the met solution
	{
#pragma omp parallel for num_threads( threads )
		for (int i = sNodes.nodeCount[depth]; i < sNodes.nodeCount[depth + 1]; i++) solution[i] = X[i - sNodes.nodeCount[depth]];
		/*
		if (!coarseToFine && depth > _minDepth)
		{
			// Explicitly compute the restriction of the met solution onto the coarser nodes
			// and down-sample the previous accumulation
			{
				UpdateConstraintsFromFiner(integrator, depth, sNodes, GetPointer(X), metConstraints + sNodes.nodeCount[depth - 1]);
				if (_constrainValues) SetPointConstraintsFromFiner(pointInfo, depth, sNodes, GetPointer(X), metConstraints + sNodes.nodeCount[depth - 1]);
				if (depth < sNodes.maxDepth - 1) DownSample(depth, sNodes, (ConstPointer(Real))metConstraints + sNodes.nodeCount[depth], metConstraints + sNodes.nodeCount[depth - 1]);
			}
		}
		*/
	}

	MemoryUsage();
	DumpOutput("\tEvaluated / Got / Solved in: %6.3f / %6.3f / %6.3f\t(%.3f MB)\n", evaluateTime, systemTime, solveTime, float(maxMemoryUsage));
	maxMemoryUsage = std::max< double >(maxMemoryUsage, _maxMemoryUsage);
	return iter;
}

template< class Real >
int Octree2< Real >::_SolveSPRGS(PointInfo& pointInfo, int depth, const typename BSplineData< 2 >::Integrator& integrator, const SortedTreeNodes& sNodes, Pointer(Real) solution, Pointer(Real) constraints, Pointer(Real) metSolutionConstraints, int iters, bool coarseToFine, bool showResidual, double* bNorm2, double* inRNorm2, double* outRNorm2, bool forceSilent)
{
	Pointer(Real) metSolution = NullPointer< Real >();
	Pointer(Real) metConstraints = NullPointer< Real >();
	if (coarseToFine) metSolution = metSolutionConstraints;	// This stores the up-sampled solution up to depth-2
	else               metConstraints = metSolutionConstraints; // This stores the down-sampled constraints up to depth

	double _maxMemoryUsage = maxMemoryUsage;
	maxMemoryUsage = 0;
	Vector< Real > X, B;
	int slices = 1 << depth;
	double systemTime = 0., solveTime = 0., updateTime = 0., evaluateTime = 0.;
	std::vector< int > offsets(slices + 1, 0);
	for (int i = sNodes.nodeCount[depth]; i < sNodes.nodeCount[depth + 1]; i++)
	{
		int d, off[3];
		sNodes.treeNodes[i]->depthAndOffset(d, off);
		offsets[off[2]]++;
	}
	for (int i = 1; i < slices; i++)  offsets[i] += offsets[i - 1];
	for (int i = slices; i >= 1; i--) offsets[i] = offsets[i - 1];
	offsets[0] = 0;

	X.Resize(sNodes.nodeCount[depth + 1] - sNodes.nodeCount[depth]);
	B.Resize(sNodes.nodeCount[depth + 1] - sNodes.nodeCount[depth]);
	if (coarseToFine)
	{
		if (depth > _minDepth)
		{
			// Up-sample the cumulative change in solution @(depth-2) into the cumulative change in solution @(depth-1)
			if (depth - 2 >= _minDepth) UpSample(depth - 1, sNodes, (ConstPointer(Real))metSolution + _sNodes.nodeCount[depth - 2], metSolution + _sNodes.nodeCount[depth - 1]);
			// Add in the change in solution @(depth-1)
#pragma omp parallel for num_threads( threads )
			for (int i = _sNodes.nodeCount[depth - 1]; i < _sNodes.nodeCount[depth]; i++) metSolution[i] += solution[i];
			// Evaluate the points @(depth) using the cumulative change in solution @(depth-1)
			/*
			if (_constrainValues)
			{
				evaluateTime = Time();
				SetPointValuesFromCoarser(pointInfo, depth, sNodes, metSolution + _sNodes.nodeCount[depth - 1]);
				evaluateTime = Time() - evaluateTime;
			}
			*/
		}
	}
	else if (depth < _sNodes.maxDepth - 1)
		for (int i = _sNodes.nodeCount[depth]; i < _sNodes.nodeCount[depth + 1]; i++) constraints[i] -= metConstraints[i];
	// Initialize with the previously computed solution
#pragma omp parallel for num_threads( threads )
	for (int i = _sNodes.nodeCount[depth]; i < _sNodes.nodeCount[depth + 1]; i++) X[i - _sNodes.nodeCount[depth]] = solution[i];
	double bNorm = 0, inRNorm = 0, outRNorm = 0;
	if (depth >= _minDepth)
	{
		int frontOffset = (showResidual || inRNorm2) ? 2 : 0;
		int backOffset = (showResidual || outRNorm2) ? 2 : 0;
		int solveSlices = std::min< int >(2 * iters - 1, slices), matrixSlices = std::max< int >(1, std::min< int >(solveSlices + frontOffset + backOffset, slices));
		std::vector< SparseMatrix< Real > > _M(matrixSlices);
		std::vector< std::vector< std::vector< int > > > __mcIndices(std::max< int >(0, solveSlices));

		int dir = coarseToFine ? -1 : 1, start = coarseToFine ? slices - 1 : 0, end = coarseToFine ? -1 : slices;
		for (int frontSlice = start - frontOffset * dir, backSlice = frontSlice - 2 * (iters - 1) * dir; backSlice != end + backOffset * dir; frontSlice += dir, backSlice += dir)
		{
			double t;
			if (frontSlice + frontOffset * dir >= 0 && frontSlice + frontOffset * dir < slices)
			{
				int s = frontSlice + frontOffset * dir, _s = s % matrixSlices;
				t = Time();
				GetSliceMatrixAndUpdateConstraints(pointInfo, _M[_s], constraints, integrator, depth, sNodes, metSolution, coarseToFine, offsets[s] + sNodes.nodeCount[depth], offsets[s + 1] + sNodes.nodeCount[depth]);
				systemTime += Time() - t;
				Pointer(TreeOctNode*) const nodes = sNodes.treeNodes + sNodes.nodeCount[depth];
				for (int i = offsets[s]; i < offsets[s + 1]; i++)
				{
					if (_boundaryType != 0 || _IsInsetSupported(nodes[i])) B[i] = constraints[nodes[i]->nodeData.nodeIndex];
					else                                                    B[i] = Real(0);
				}
				if (showResidual || inRNorm2)
#pragma omp parallel for num_threads( threads ) reduction( + : bNorm , inRNorm )
					for (int j = 0; j < _M[_s].rows; j++)
					{
						Real temp = Real(0);
						ConstPointer(MatrixEntry< Real >) start = _M[_s][j];
						ConstPointer(MatrixEntry< Real >) end = start + _M[_s].rowSizes[j];
						ConstPointer(MatrixEntry< Real >) e;
						for (e = start; e != end; e++) temp += X[e->N] * e->Value;
						Real b = B[j + offsets[s]];
						bNorm += b * b;
						inRNorm += (temp - b) * (temp - b);
					}
				else if (bNorm2)
#pragma omp parallel for num_threads( threads ) reduction( + : bNorm )
					for (int j = 0; j < _M[_s].rows; j++)
					{
						Real b = B[j + offsets[s]];
						bNorm += b * b;
					}
			}
			t = Time();
			if (iters && frontSlice >= 0 && frontSlice < slices)
			{
				int s = frontSlice, _s = s % matrixSlices, __s = s % solveSlices;
				for (int i = 0; i<int(__mcIndices[__s].size()); i++) __mcIndices[__s][i].clear();
				_setMultiColorIndices(sNodes.nodeCount[depth] + offsets[s], sNodes.nodeCount[depth] + offsets[s + 1], __mcIndices[__s]);
			}
			for (int slice = frontSlice; slice * dir >= backSlice * dir; slice -= 2 * dir)
				if (slice >= 0 && slice < slices)
				{
					int s = slice, _s = s % matrixSlices, __s = s % solveSlices;
					SparseMatrix< Real >::SolveGS(__mcIndices[__s], _M[_s], B, X, !coarseToFine, threads, offsets[s]);
				}
			solveTime += Time() - t;
			if ((showResidual || outRNorm2) && backSlice - backOffset * dir >= 0 && backSlice - backOffset * dir < slices)
			{
				int s = backSlice - backOffset * dir, _s = s % matrixSlices;
#pragma omp parallel for num_threads( threads ) reduction( + : outRNorm )
				for (int j = 0; j < _M[_s].rows; j++)
				{
					Real temp = Real(0);
					ConstPointer(MatrixEntry< Real >) start = _M[_s][j];
					ConstPointer(MatrixEntry< Real >) end = start + _M[_s].rowSizes[j];
					ConstPointer(MatrixEntry< Real >) e;
					for (e = start; e != end; e++) temp += X[e->N] * e->Value;
					Real b = B[j + offsets[s]];
					outRNorm += (temp - b) * (temp - b);
				}
			}
		}
	}

	if (bNorm2) bNorm2[depth] = bNorm;
	if (inRNorm2) inRNorm2[depth] = inRNorm;
	if (outRNorm2) outRNorm2[depth] = outRNorm;
	if (showResidual && iters)
	{
		for (int i = 0; i < depth; i++) printf("  ");
		printf("GS: %.4e -> %.4e -> %.4e (%.2e) [%d]\n", sqrt(bNorm), sqrt(inRNorm), sqrt(outRNorm), sqrt(outRNorm / bNorm), iters);
	}

	// Copy the old solution into the buffer, write in the new solution, compute the change, and update the met constraints
#pragma omp parallel for num_threads( threads )
	for (int i = sNodes.nodeCount[depth]; i < sNodes.nodeCount[depth + 1]; i++) solution[i] = X[i - sNodes.nodeCount[depth]];
	if (!coarseToFine && depth > _minDepth)
	{
		// Explicitly compute the restriction of the met solution onto the coarser nodes
		// and down-sample the previous accumulation
		{
			UpdateConstraintsFromFiner(integrator, depth, sNodes, GetPointer(X), metConstraints + sNodes.nodeCount[depth - 1]);
			if (_constrainValues) SetPointConstraintsFromFiner(pointInfo, depth, sNodes, GetPointer(X), metConstraints + sNodes.nodeCount[depth - 1]);
			if (depth < sNodes.maxDepth - 1) DownSample(depth, sNodes, (ConstPointer(Real))metConstraints + sNodes.nodeCount[depth], metConstraints + sNodes.nodeCount[depth - 1]);
		}
	}

	MemoryUsage();
	if (!forceSilent) DumpOutput("\tEvaluated / Got / Solved in: %6.3f / %6.3f / %6.3f\t(%.3f MB)\n", evaluateTime, systemTime, solveTime, float(maxMemoryUsage));
	maxMemoryUsage = std::max< double >(maxMemoryUsage, _maxMemoryUsage);

	return iters;
}

template< class Real >
int Octree2< Real >::HasNormals(TreeOctNode* node, const SmoothCoeff& smoothcoeff)
{
	int idx = smoothcoeff.scoeffIndex(node);
	if (idx >= 0)
	{
		
		if (smoothcoeff.weight[idx] != 0 
			&& (smoothcoeff.boundary[idx][0] == false 
				|| smoothcoeff.boundary[idx][1] == false 
				   || smoothcoeff.boundary[idx][2] == false))
		{
			return 1;
		}
		/*
		if (smoothcoeff.weight[idx] != 0)
		{
			return 1;
		}
		*/
	}
	if (node->children) for (int i = 0; i < Cube::CORNERS; i++) if (HasNormals(&node->children[i], smoothcoeff)) return 1;
	return 0;
}
template< class Real >
std::vector < std::vector< MatrixEntry <Real> >>  Octree2< Real >::SetLaplacianConstraints2(const SmoothCoeff& smoothcoeff)
{
	typename BSplineData< 2 >::Integrator integrator;
	_fData.setIntegrator(integrator, _boundaryType == 0);
	int maxDepth = _sNodes.maxDepth - 1;
	Point3D< Real > zeroPoint;
	zeroPoint[0] = zeroPoint[1] = zeroPoint[2] = 0;
	Real thre = 1e-6;
	std::vector < std::vector< MatrixEntry <Real> > > Q(_sNodes.nodeCount[_sNodes.maxDepth], std::vector< MatrixEntry <Real> >());
	std::vector < std::unordered_map<int, int>> Q_index(_sNodes.nodeCount[_sNodes.maxDepth]); //Matrix index to vector index 
	std::vector < std::vector< MatrixEntry <Real> > > _Q(_sNodes.nodeCount[maxDepth], std::vector< MatrixEntry <Real> >());
	std::vector < std::unordered_map<int, int>> _Q_index(_sNodes.nodeCount[maxDepth]); 
	
	MemoryUsage();
	//Constructing the matrix Q (named B in the paper)
	for (int d = maxDepth; d >= (_boundaryType == 0 ? 2 : 0); d--)
	{
		int offset = d > 0 ? _sNodes.treeNodes[_sNodes.nodeCount[d - 1]]->nodeData.nodeIndex : 0;
		Stencil< Point3D< double >, 5 > stencil, stencils[2][2][2];
		
		SetDivergenceStencil(d, integrator, stencil, false);
		SetDivergenceStencils(d, integrator, stencils, true);

		std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
		for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(_fData.depth);
#pragma omp parallel for num_threads( threads )
		for (int i = _sNodes.nodeCount[d]; i < _sNodes.nodeCount[d + 1]; i++)
		{
			typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
			TreeOctNode* node = _sNodes.treeNodes[i];
			int startX = 0, endX = 5, startY = 0, endY = 5, startZ = 0, endZ = 5;
			int depth = node->depth();
			typename TreeOctNode::Neighbors5 neighbors5;
			neighborKey.getNeighbors(node, neighbors5);

			bool isInterior, isInterior2;
			{
				int d, off[3];
				node->depthAndOffset(d, off);
				int o = _boundaryType == 0 ? (1 << (d - 2)) : 0;
				int mn = 2 + o, mx = (1 << d) - 2 - o;
				isInterior = (off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx);
				mn += 2, mx -= 2;
				isInterior2 = (off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx);
			}
			int cx, cy, cz;
			if (d)
			{
				int c = int(node - node->parent->children);
				Cube::FactorCornerIndex(c, cx, cy, cz);
			}
			else cx = cy = cz = 0;
			Stencil< Point3D< double >, 5 >& _stencil = stencils[cx][cy][cz];
			int d, off[3];
			node->depthAndOffset(d, off);
			// Set constraints from current depth
			// Gather the constraints from the vector-field at _node into the constraint stored with node
			{
				if (isInterior)
				{
					for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
					{
						const TreeOctNode* _node = neighbors5.neighbors[x][y][z];
						if (_node)
						{
							int _idx = smoothcoeff.scoeffIndex(_node);
							if (_idx >= 0)
							{
								std::vector < MatrixEntry <Real>> Q_new_line;
								const std::vector<MatrixEntry <Real>>& _coeff = smoothcoeff.smoothcoeff[_idx];
								int origin_size = _coeff.size();
								Q_new_line.resize(3 * origin_size);
								//printf("%d\n", smoothcoeff.smoothcoeff[_idx].size());
								//3n dim
								for (int j = 0; j < origin_size; j++)
								{
									Real value0 = Real(stencil.values[x][y][z].coords[0] * _coeff[j].Value);
									Real value1 = Real(stencil.values[x][y][z].coords[1] * _coeff[j].Value);
									Real value2 = Real(stencil.values[x][y][z].coords[2] * _coeff[j].Value);
									if (value0 != 0 && smoothcoeff.boundary[_idx][0] == false) //不是边界
									{
										Q_new_line[3 * j] = MatrixEntry< Real >(_coeff[j].N * 3, value0);
									}
									else
									{
										Q_new_line[3 * j] = MatrixEntry< Real >(_coeff[j].N * 3, 0);
									}
									if (value1 != 0 && smoothcoeff.boundary[_idx][1] == false)
									{
										Q_new_line[3 * j + 1] = MatrixEntry< Real >(_coeff[j].N * 3 + 1, value1);
									}
									else
									{
										Q_new_line[3 * j + 1] = MatrixEntry< Real >(_coeff[j].N * 3 + 1, 0);
									}
									if (value2 != 0 && smoothcoeff.boundary[_idx][2] == false)
									{
										Q_new_line[3 * j + 2] = MatrixEntry< Real >(_coeff[j].N * 3 + 2, value2);
									}
									else
									{
										Q_new_line[3 * j + 2] = MatrixEntry< Real >(_coeff[j].N * 3 + 2, 0);
									}
									//Point3D<Real> n = all_normals[smoothcoeff.smoothcoeff[_idx][j].N];
								}
								//printf("%d\n", Q_new_line.size());
								if (Q_new_line.size() > 0)
								{
									//printf("%d\n", Q[node->nodeData.nodeIndex].size());
									AddSparseVector2(Q[node->nodeData.nodeIndex], Q_index[node->nodeData.nodeIndex], Q_new_line, thre);
								}
							}
						}
					}
				}
				else
				{
					for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
					{
						const TreeOctNode* _node = neighbors5.neighbors[x][y][z];
						if (_node)
						{
							int _idx = smoothcoeff.scoeffIndex(_node);
							if (_idx >= 0)
							{
								int _d, _off[3];
								_node->depthAndOffset(_d, _off);
								std::vector < MatrixEntry <Real>> Q_new_line;
								Point3D<double> &div = GetDivergence2(integrator, d, off, _off, false);
								const std::vector<MatrixEntry <Real>>& _coeff = smoothcoeff.smoothcoeff[_idx];
								int origin_size = _coeff.size();
								Q_new_line.resize(3 * origin_size);
								for (int j = 0; j < origin_size; j++)
								{
									Real value0 = Real(div.coords[0] * _coeff[j].Value);
									Real value1 = Real(div.coords[1] * _coeff[j].Value);
									Real value2 = Real(div.coords[2] * _coeff[j].Value);
									if (value0 != 0 && smoothcoeff.boundary[_idx][0] == false) //不是边界
									{
										Q_new_line[3 * j] = MatrixEntry< Real >(_coeff[j].N * 3, value0);
									}
									else
									{
										Q_new_line[3 * j] = MatrixEntry< Real >(_coeff[j].N * 3, 0);
									}
									if (value1 != 0 && smoothcoeff.boundary[_idx][1] == false)
									{
										Q_new_line[3 * j + 1] = MatrixEntry< Real >(_coeff[j].N * 3 + 1, value1);
									}
									else
									{
										Q_new_line[3 * j + 1] = MatrixEntry< Real >(_coeff[j].N * 3 + 1, 0);
									}
									if (value2 != 0 && smoothcoeff.boundary[_idx][2] == false)
									{
										Q_new_line[3 * j + 2] = MatrixEntry< Real >(_coeff[j].N * 3 + 2, value2);
									}
									else
									{
										Q_new_line[3 * j + 2] = MatrixEntry< Real >(_coeff[j].N * 3 + 2, 0);
									}
									//Point3D<Real> n = all_normals[smoothcoeff.smoothcoeff[_idx][j].N];
								}
								//printf("%d\n", Q_new_line.size());
								if (Q_new_line.size() > 0)
								{
									//printf("%d\n", Q[node->nodeData.nodeIndex].size());
									AddSparseVector2(Q[node->nodeData.nodeIndex], Q_index[node->nodeData.nodeIndex], Q_new_line, thre);
								}
							}
						}
					}
				}
				UpdateCoarserSupportBounds(neighbors5.neighbors[2][2][2], startX, endX, startY, endY, startZ, endZ);
			}
			int idx = smoothcoeff.scoeffIndex(node);
			if (idx < 0) continue;
			const std::vector< MatrixEntry <Real> >& coeff = smoothcoeff.smoothcoeff[idx];
			const std::vector<bool>& bound = smoothcoeff.boundary[idx];
			if (smoothcoeff.weight[idx] == 0 || (bound[0] == bound[1] == bound[2] == true)) continue;

			// Set the constraints for the parents
			if (depth > _minDepth)
			{
				neighborKey.getNeighbors(node->parent, neighbors5);

				for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
				{
					if (neighbors5.neighbors[x][y][z])
					{
						TreeOctNode* _node = neighbors5.neighbors[x][y][z];
						std::vector < MatrixEntry <Real>> Q_new_line;
						//Real c;
						if (isInterior2)
						{
							Point3D< double >& div = _stencil.values[x][y][z];
							int origin_size = coeff.size();
							Q_new_line.resize(3 * origin_size);
							for (int j = 0; j < origin_size; j++)
							{
								Real value0 = Real(div.coords[0] * coeff[j].Value);
								Real value1 = Real(div.coords[1] * coeff[j].Value);
								Real value2 = Real(div.coords[2] * coeff[j].Value);
								if (value0 != 0 && bound[0] == false)
								{
									Q_new_line[3 * j] = MatrixEntry< Real >(coeff[j].N * 3, value0);
								}
								else
								{
									Q_new_line[3 * j] = MatrixEntry< Real >(coeff[j].N * 3, 0);
								}
								if (value1 != 0 && bound[1] == false)
								{
									Q_new_line[3 * j + 1] = MatrixEntry< Real >(coeff[j].N * 3 + 1, value1);
								}
								else
								{
									Q_new_line[3 * j + 1] = MatrixEntry< Real >(coeff[j].N * 3 + 1, 0);
								}
								if (value2 != 0 && bound[2] == false)
								{
									Q_new_line[3 * j + 2] = MatrixEntry< Real >(coeff[j].N * 3 + 2, value2);
								}
								else
								{
									Q_new_line[3 * j + 2] = MatrixEntry< Real >(coeff[j].N * 3 + 2, 0);
								}
							}
							//c = Real(div[0] * normal[0] + div[1] * normal[1] + div[2] * normal[2]);
						}
						else
						{
							int _d, _off[3];
							_node->depthAndOffset(_d, _off);
							Point3D<double>& div = GetDivergence1(integrator, d, off, _off, true);
							int origin_size = coeff.size();
							Q_new_line.resize(3 * origin_size);
							for (int j = 0; j < origin_size; j++)
							{
								Real value0 = Real(div[0] * coeff[j].Value);
								Real value1 = Real(div[1] * coeff[j].Value);
								Real value2 = Real(div[2] * coeff[j].Value);
								if (value0 != 0 && bound[0] == false)
								{
									Q_new_line[3 * j] = MatrixEntry< Real >(coeff[j].N * 3, value0);
								}
								else
								{
									Q_new_line[3 * j] = MatrixEntry< Real >(coeff[j].N * 3, 0);
								}
								if (value1 != 0 && bound[1] == false)
								{
									Q_new_line[3 * j + 1] = MatrixEntry< Real >(coeff[j].N * 3 + 1, value1);
								}
								else
								{
									Q_new_line[3 * j + 1] = MatrixEntry< Real >(coeff[j].N * 3 + 1, 0);
								}
								if (value2 != 0 && bound[2] == false)
								{
									Q_new_line[3 * j + 2] = MatrixEntry< Real >(coeff[j].N * 3 + 2, value2);
								}
								else
								{
									Q_new_line[3 * j + 2] = MatrixEntry< Real >(coeff[j].N * 3 + 2, 0);
								}
							}
						}
#pragma omp critical
						{
							if (Q_new_line.size() > 0)
							{
								//printf("%d\n", _Q[_node->nodeData.nodeIndex].size());
								AddSparseVector2(_Q[_node->nodeData.nodeIndex], _Q_index[_node->nodeData.nodeIndex], Q_new_line, thre);
							}
						}
					}
				}
			}
		}

		//printf("%d\n", d);
		MemoryUsage();
	}

	

	// Fine-to-coarse down-sampling of constraints
	//for (int d = maxDepth - 1; d >= (_boundaryType == 0 ? 2 : 0); d--) DownSample(d, _sNodes, (ConstPointer(Real))_constraints + _sNodes.nodeCount[d], _constraints + _sNodes.nodeCount[d - 1]);
	for (int d = maxDepth - 1; d >= (_boundaryType == 0 ? 2 : 0); d--) DownSample3(d, _sNodes, _sNodes.nodeCount[d], _sNodes.nodeCount[d - 1], _Q, _Q_index, thre);

	// Add the accumulated constraints from all finer depths
#pragma omp parallel for num_threads( threads )
	for (int i = 0; i < _sNodes.nodeCount[maxDepth]; i++) AddSparseVector2(Q[i], Q_index[i], _Q[i], thre);

	//FreePointer(_constraints);


	
	//relation between sum normals to each sample point
	std::vector< std::vector < MatrixEntry <Real>> > coefficients1(_sNodes.nodeCount[maxDepth], std::vector < MatrixEntry <Real>>());
	std::vector < std::unordered_map<int, int>> coefficients_index1(_sNodes.nodeCount[maxDepth]);
	std::vector< std::vector < MatrixEntry <Real>> > coefficients2(_sNodes.nodeCount[maxDepth], std::vector < MatrixEntry <Real>>());
	std::vector < std::unordered_map<int, int>> coefficients_index2(_sNodes.nodeCount[maxDepth]);
	std::vector< std::vector < MatrixEntry <Real>> > coefficients3(_sNodes.nodeCount[maxDepth], std::vector < MatrixEntry <Real>>());
	std::vector < std::unordered_map<int, int>> coefficients_index3(_sNodes.nodeCount[maxDepth]);
	std::vector< std::vector < bool> > boundarys(_sNodes.nodeCount[maxDepth], std::vector < bool>(3, false));
	for (int d = maxDepth - 1; d >= 0; d--)
	{
#pragma omp parallel for num_threads( threads )
		for (int i = _sNodes.nodeCount[d]; i < _sNodes.nodeCount[d + 1]; i++)
		{
			int idx = smoothcoeff.scoeffIndex(_sNodes.treeNodes[i]);
			if (idx < 0) continue;
			int origin_size = smoothcoeff.smoothcoeff[idx].size();
			for (int j = 0; j < 3; j++)
			{
				boundarys[i][j] = smoothcoeff.boundary[idx][j];
			}
			for (int j = 0; j < origin_size; j++)
			{
				if (boundarys[i][0] == false)
				{
					std::vector < MatrixEntry <Real>>& coeff = coefficients1[i];
					coeff.push_back(smoothcoeff.smoothcoeff[idx][j]);
					std::unordered_map<int, int>& coefficients_map = coefficients_index1[i];
					std::pair<int, int> pair((int)smoothcoeff.smoothcoeff[idx][j].N, j);
					coefficients_map.insert(pair);
				}
				if (boundarys[i][1] == false)
				{
					std::vector < MatrixEntry <Real>>& coeff = coefficients2[i];
					coeff.push_back(smoothcoeff.smoothcoeff[idx][j]);
					std::unordered_map<int, int>& coefficients_map = coefficients_index2[i];
					std::pair<int, int> pair((int)smoothcoeff.smoothcoeff[idx][j].N, j);
					coefficients_map.insert(pair);
				}
				if (boundarys[i][2] == false)
				{
					std::vector < MatrixEntry <Real>>& coeff = coefficients3[i];
					coeff.push_back(smoothcoeff.smoothcoeff[idx][j]);
					std::unordered_map<int, int>& coefficients_map = coefficients_index3[i];
					std::pair<int, int> pair((int)smoothcoeff.smoothcoeff[idx][j].N, j);
					coefficients_map.insert(pair);
				}
			}
		}
	}
	
	// Coarse-to-fine up-sampling of coefficients
	for (int d = (_boundaryType == 0 ? 2 : 0); d < maxDepth; d++) UpSample3(d, _sNodes, _sNodes.nodeCount[d - 1], _sNodes.nodeCount[d], coefficients1, coefficients_index1, thre);
	for (int d = (_boundaryType == 0 ? 2 : 0); d < maxDepth; d++) UpSample3(d, _sNodes, _sNodes.nodeCount[d - 1], _sNodes.nodeCount[d], coefficients2, coefficients_index2, thre);
	for (int d = (_boundaryType == 0 ? 2 : 0); d < maxDepth; d++) UpSample3(d, _sNodes, _sNodes.nodeCount[d - 1], _sNodes.nodeCount[d], coefficients3, coefficients_index3, thre);


	// Compute the contribution from all coarser depths
	for (int d = 0; d <= maxDepth; d++)
	{
		size_t start = _sNodes.nodeCount[d], end = _sNodes.nodeCount[d + 1], range = end - start;
		Stencil< Point3D< double >, 5 > stencils[2][2][2];
		SetDivergenceStencils(d, integrator, stencils, false);
		std::vector< typename TreeOctNode::NeighborKey3 > neighborKeys(std::max< int >(1, threads));
		for (int i = 0; i < neighborKeys.size(); i++) neighborKeys[i].set(maxDepth);
#pragma omp parallel for num_threads( threads )
		for (int i = _sNodes.nodeCount[d]; i < _sNodes.nodeCount[d + 1]; i++)
		{
			typename TreeOctNode::NeighborKey3& neighborKey = neighborKeys[omp_get_thread_num()];
			TreeOctNode* node = _sNodes.treeNodes[i];
			int depth = node->depth();
			if (!depth) continue;
			int startX = 0, endX = 5, startY = 0, endY = 5, startZ = 0, endZ = 5;
			UpdateCoarserSupportBounds(node, startX, endX, startY, endY, startZ, endZ);
			typename TreeOctNode::Neighbors5 neighbors5;
			neighborKey.getNeighbors(node->parent, neighbors5);

			bool isInterior;
			{
				int d, off[3];
				node->depthAndOffset(d, off);
				int o = _boundaryType == 0 ? (1 << (d - 2)) : 0;
				int mn = 4 + o, mx = (1 << d) - 4 - o;
				isInterior = (off[0] >= mn && off[0] < mx && off[1] >= mn && off[1] < mx && off[2] >= mn && off[2] < mx);
			}
			int cx, cy, cz;
			if (d)
			{
				int c = int(node - node->parent->children);
				Cube::FactorCornerIndex(c, cx, cy, cz);
			}
			else cx = cy = cz = 0;
			Stencil< Point3D< double >, 5 >& _stencil = stencils[cx][cy][cz];

			//Real constraint = Real(0);
			int d, off[3];
			node->depthAndOffset(d, off);
			std::vector < MatrixEntry <Real>> constraint_one_line(0);
			std::unordered_map<int, int> constraint_one_line_index;
			for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
			{
				if (neighbors5.neighbors[x][y][z])
				{
					TreeOctNode* _node = neighbors5.neighbors[x][y][z];
					int _i = _node->nodeData.nodeIndex;
					Point3D< double > div;
					if (isInterior)
					{
						div = _stencil.values[x][y][z];
						//Point3D< Real >& normal = coefficients[_i];
						//constraint += Real(div[0] * normal[0] + div[1] * normal[1] + div[2] * normal[2]);
					}
					else
					{
						int _d, _off[3];
						_node->depthAndOffset(_d, _off);
						div = GetDivergence2(integrator, d, off, _off, true);
					}
					std::vector< MatrixEntry <Real> >& coeff1 = coefficients1[_i];
					std::vector< MatrixEntry <Real> >& coeff2 = coefficients2[_i];
					std::vector< MatrixEntry <Real> >& coeff3 = coefficients3[_i];
					std::vector < MatrixEntry <Real>> Q_new_line;
					int origin_size1 = coeff1.size();
					int origin_size2 = coeff2.size();
					int origin_size3 = coeff3.size();
					Q_new_line.resize(origin_size1 + origin_size2 + origin_size3);
					for (int j = 0; j < origin_size1; j++)
					{
						Real value0 = Real(div[0] * coeff1[j].Value);
						if (value0 != 0)
						{
							Q_new_line[j] = MatrixEntry< Real >(coeff1[j].N * 3, value0);
						}
					}
					for (int j = 0; j < origin_size2; j++)
					{
						Real value1 = Real(div[1] * coeff2[j].Value);
						if (value1 != 0)
						{
							Q_new_line[j + origin_size1] = MatrixEntry< Real >(coeff2[j].N * 3 + 1, value1);
						}
					}
					for (int j = 0; j < origin_size3; j++)
					{
						Real value2 = Real(div[2] * coeff3[j].Value);
						if (value2 != 0)
						{
							Q_new_line[j + origin_size1 + origin_size2] = MatrixEntry< Real >(coeff3[j].N * 3 + 2, value2);
						}
					}
					if (Q_new_line.size() > 0)
					{
						//printf("%d\n", Q[node->nodeData.nodeIndex].size());
						AddSparseVector2(constraint_one_line, constraint_one_line_index, Q_new_line, thre);
					}
				}
			}
			AddSparseVector2(Q[node->nodeData.nodeIndex], Q_index[node->nodeData.nodeIndex], constraint_one_line, thre);
		}
	}
	MemoryUsage();
	return Q;
}
template< class Real >
void Octree2< Real >::refineBoundary(std::vector< int >* map) { _sNodes.set(tree, map); }



template< class Real >
Real Octree2< Real >::getCenterValue(const typename TreeOctNode::ConstNeighborKey3& neighborKey, const TreeOctNode* node, ConstPointer(Real) solution, ConstPointer(Real) metSolution, const typename BSplineData< 2 >::template CenterEvaluator< 1 >& evaluator, const Stencil< double, 3 >& stencil, const Stencil< double, 3 >& pStencil, bool isInterior) const
{
	if (node->children) fprintf(stderr, "[WARNING] getCenterValue assumes leaf node\n");
	Real value = 0;

	int d, off[3];
	node->depthAndOffset(d, off);

	if (isInterior)
	{
		for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++)
		{
			const TreeOctNode* n = neighborKey.neighbors[d].neighbors[i][j][k];
			if (n) value += solution[n->nodeData.nodeIndex] * Real(stencil.values[i][j][k]);
		}
		if (d > _minDepth)
			for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++)
			{
				const TreeOctNode* n = neighborKey.neighbors[d - 1].neighbors[i][j][k];
				if (n) value += metSolution[n->nodeData.nodeIndex] * Real(pStencil.values[i][j][k]);
			}
	}
	else
	{
		for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++)
		{
			const TreeOctNode* n = neighborKey.neighbors[d].neighbors[i][j][k];
			if (n)
			{
				int _d, _off[3];
				n->depthAndOffset(_d, _off);
				value +=
					solution[n->nodeData.nodeIndex] * Real(
						evaluator.value(d, off[0], _off[0], false, false) * evaluator.value(d, off[1], _off[1], false, false) * evaluator.value(d, off[1], _off[1], false, false));
			}
		}
		if (d > _minDepth)
			for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++)
			{
				const TreeOctNode* n = neighborKey.neighbors[d - 1].neighbors[i][j][k];
				if (n)
				{
					int _d, _off[3];
					n->depthAndOffset(_d, _off);
					value +=
						solution[n->nodeData.nodeIndex] * Real(
							evaluator.value(d, off[0], _off[0], false, false) * evaluator.value(d, off[1], _off[1], false, false) * evaluator.value(d, off[1], _off[1], false, false));
				}
			}
	}
	return value;
}
template< class Real >
Real Octree2< Real >::getCornerValue(const typename TreeOctNode::ConstNeighborKey3& neighborKey, const TreeOctNode* node, int corner, ConstPointer(Real) solution, ConstPointer(Real) metSolution, const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator, const Stencil< double, 3 >& stencil, const Stencil< double, 3 > stencils[8], bool isInterior) const
{
	double value = 0;
	if (_boundaryType == -1) value = -0.5;
	int d, off[3];
	node->depthAndOffset(d, off);

	int cx, cy, cz;
	int startX = 0, endX = 3, startY = 0, endY = 3, startZ = 0, endZ = 3;
	Cube::FactorCornerIndex(corner, cx, cy, cz);
	{
		typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.neighbors[d];
		if (cx == 0) endX = 2;
		else      startX = 1;
		if (cy == 0) endY = 2;
		else      startY = 1;
		if (cz == 0) endZ = 2;
		else      startZ = 1;
		if (isInterior)
			for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if (_node) value += solution[_node->nodeData.nodeIndex] * stencil.values[x][y][z];
			}
		else
			for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if (_node)
				{
					int _d, _off[3];
					_node->depthAndOffset(_d, _off);
					value += solution[_node->nodeData.nodeIndex] * evaluator.value(d, off[0], cx, _off[0], false, false) * evaluator.value(d, off[1], cy, _off[1], false, false) * evaluator.value(d, off[2], cz, _off[2], false, false);
				}
			}
	}
	if (d > _minDepth)
	{
		int _corner = int(node - node->parent->children);
		int _cx, _cy, _cz;
		Cube::FactorCornerIndex(_corner, _cx, _cy, _cz);
		if (cx != _cx) startX = 0, endX = 3;
		if (cy != _cy) startY = 0, endY = 3;
		if (cz != _cz) startZ = 0, endZ = 3;
		typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.neighbors[d - 1];
		if (isInterior)
			for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if (_node) value += metSolution[_node->nodeData.nodeIndex] * stencils[_corner].values[x][y][z];
			}
		else
			for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if (_node)
				{
					int _d, _off[3];
					_node->depthAndOffset(_d, _off);
					value += metSolution[_node->nodeData.nodeIndex] * evaluator.value(d, off[0], cx, _off[0], false, true) * evaluator.value(d, off[1], cy, _off[1], false, true) * evaluator.value(d, off[2], cz, _off[2], false, true);
				}
			}
	}
	return Real(value);
}
template< class Real >
std::pair< Real, Point3D< Real > > Octree2< Real >::getCornerValueAndNormal(const typename TreeOctNode::ConstNeighborKey3& neighborKey, const TreeOctNode* node, int corner, ConstPointer(Real) solution, ConstPointer(Real) metSolution, const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator, const Stencil< double, 3 >& vStencil, const Stencil< double, 3 > vStencils[8], const Stencil< Point3D< double >, 3 >& nStencil, const Stencil< Point3D< double >, 3 > nStencils[8], bool isInterior) const
{
	double value = 0;
	Point3D< double > normal;
	if (_boundaryType == -1) value = -0.5;
	int d, off[3];
	node->depthAndOffset(d, off);

	int cx, cy, cz;
	int startX = 0, endX = 3, startY = 0, endY = 3, startZ = 0, endZ = 3;
	Cube::FactorCornerIndex(corner, cx, cy, cz);
	{
		typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.neighbors[d];
		if (cx == 0) endX = 2;
		else      startX = 1;
		if (cy == 0) endY = 2;
		else      startY = 1;
		if (cz == 0) endZ = 2;
		else      startZ = 1;
		if (isInterior)
			for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if (_node) value += solution[_node->nodeData.nodeIndex] * vStencil.values[x][y][z], normal += nStencil.values[x][y][z] * solution[_node->nodeData.nodeIndex];
			}
		else
			for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if (_node)
				{
					int _d, _off[3];
					_node->depthAndOffset(_d, _off);
					double v[] = { evaluator.value(d , off[0] , cx , _off[0] , false , false) , evaluator.value(d , off[1] , cy , _off[1] , false , false) , evaluator.value(d , off[2] , cz , _off[2] , false , false) };
					double dv[] = { evaluator.value(d , off[0] , cx , _off[0] , true  , false) , evaluator.value(d , off[1] , cy , _off[1] , true  , false) , evaluator.value(d , off[2] , cz , _off[2] , true  , false) };
					value += solution[_node->nodeData.nodeIndex] * evaluator.value(d, off[0], cx, _off[0], false, false) * evaluator.value(d, off[1], cy, _off[1], false, false) * evaluator.value(d, off[2], cz, _off[2], false, false);
					normal += Point3D< double >(dv[0] * v[1] * v[2], v[0] * dv[1] * v[2], v[0] * v[1] * dv[2]) * solution[_node->nodeData.nodeIndex];
				}
			}
	}
	if (d > _minDepth)
	{
		int _corner = int(node - node->parent->children);
		int _cx, _cy, _cz;
		Cube::FactorCornerIndex(_corner, _cx, _cy, _cz);
		if (cx != _cx) startX = 0, endX = 3;
		if (cy != _cy) startY = 0, endY = 3;
		if (cz != _cz) startZ = 0, endZ = 3;
		typename TreeOctNode::ConstNeighbors3& neighbors = neighborKey.neighbors[d - 1];
		if (isInterior)
			for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if (_node) value += metSolution[_node->nodeData.nodeIndex] * vStencils[_corner].values[x][y][z], normal += nStencils[_corner].values[x][y][z] * metSolution[_node->nodeData.nodeIndex];
			}
		else
			for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if (_node)
				{
					int _d, _off[3];
					_node->depthAndOffset(_d, _off);
					double v[] = { evaluator.value(d , off[0] , cx , _off[0] , false , true) , evaluator.value(d , off[1] , cy , _off[1] , false , true) , evaluator.value(d , off[2] , cz , _off[2] , false , true) };
					double dv[] = { evaluator.value(d , off[0] , cx , _off[0] , true  , true) , evaluator.value(d , off[1] , cy , _off[1] , true  , true) , evaluator.value(d , off[2] , cz , _off[2] , true  , true) };
					value += metSolution[_node->nodeData.nodeIndex] * evaluator.value(d, off[0], cx, _off[0], false, true) * evaluator.value(d, off[1], cy, _off[1], false, true) * evaluator.value(d, off[2], cz, _off[2], false, true);
					normal += Point3D< double >(dv[0] * v[1] * v[2], v[0] * dv[1] * v[2], v[0] * v[1] * dv[2]) * metSolution[_node->nodeData.nodeIndex];
				}
			}
	}
	return std::pair< Real, Point3D< Real > >(Real(value), Point3D< Real >(normal));
}
template< class Real >
Point3D< Real > Octree2< Real >::getCornerNormal(const typename TreeOctNode::ConstNeighbors5& neighbors5, const typename TreeOctNode::ConstNeighbors5& pNeighbors5, const TreeOctNode* node, int corner, ConstPointer(Real) solution, ConstPointer(Real) metSolution, const typename BSplineData< 2 >::template CornerEvaluator< 2 >& evaluator, const Stencil< Point3D< double >, 5 >& nStencil, const Stencil< Point3D< double >, 5 > nStencils[8], bool isInterior) const
{
	Point3D< double > normal;
	normal[0] = normal[1] = normal[2] = 0.;

	int d, off[3];
	node->depthAndOffset(d, off);

	int cx, cy, cz;
	int startX = 0, endX = 5, startY = 0, endY = 5, startZ = 0, endZ = 5;
	Cube::FactorCornerIndex(corner, cx, cy, cz);
	{
		if (cx == 0) endX = 4;
		else      startX = 1;
		if (cy == 0) endY = 4;
		else      startY = 1;
		if (cz == 0) endZ = 4;
		else      startZ = 1;
		if (isInterior)
			for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
			{
				const TreeOctNode* _node = neighbors5.neighbors[x][y][z];
				if (_node) normal += nStencil.values[x][y][z] * solution[_node->nodeData.nodeIndex];
			}
		else
			for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
			{
				const TreeOctNode* _node = neighbors5.neighbors[x][y][z];
				if (_node)
				{
					int _d, _off[3];
					_node->depthAndOffset(_d, _off);
					double v[] = { evaluator.value(d , off[0] , cx , _off[0] , false , false) , evaluator.value(d , off[1] , cy , _off[1] , false , false) , evaluator.value(d , off[2] , cz , _off[2] , false , false) };
					double dv[] = { evaluator.value(d , off[0] , cx , _off[0] , true  , false) , evaluator.value(d , off[1] , cy , _off[1] , true  , false) , evaluator.value(d , off[2] , cz , _off[2] , true  , false) };
					normal += Point3D< double >(dv[0] * v[1] * v[2], v[0] * dv[1] * v[2], v[0] * v[1] * dv[2]) * solution[_node->nodeData.nodeIndex];
				}
			}
	}
	if (d > _minDepth)
	{
		int _cx, _cy, _cz, _corner = int(node - node->parent->children);
		Cube::FactorCornerIndex(_corner, _cx, _cy, _cz);
		if (cx != _cx) startX = 0, endX = 5;
		if (cy != _cy) startY = 0, endY = 5;
		if (cz != _cz) startZ = 0, endZ = 5;
		if (isInterior)
			for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
			{
				const TreeOctNode* _node = pNeighbors5.neighbors[x][y][z];
				if (_node) normal += nStencils[_corner].values[x][y][z] * metSolution[_node->nodeData.nodeIndex];
			}
		else
			for (int x = startX; x < endX; x++) for (int y = startY; y < endY; y++) for (int z = startZ; z < endZ; z++)
			{
				const TreeOctNode* _node = pNeighbors5.neighbors[x][y][z];
				if (_node)
				{
					int _d, _off[3];
					_node->depthAndOffset(_d, _off);
					double v[] = { evaluator.value(d , off[0] , cx , _off[0] , false , true) , evaluator.value(d , off[1] , cy , _off[1] , false , true) , evaluator.value(d , off[2] , cz , _off[2] , false , true) };
					double dv[] = { evaluator.value(d , off[0] , cx , _off[0] , true  , true) , evaluator.value(d , off[1] , cy , _off[1] , true  , true) , evaluator.value(d , off[2] , cz , _off[2] , true  , true) };
					normal += Point3D< double >(dv[0] * v[1] * v[2], v[0] * dv[1] * v[2], v[0] * v[1] * dv[2]) * metSolution[_node->nodeData.nodeIndex];
				}
			}
	}
	return Point3D< Real >(Real(normal[0]), Real(normal[1]), Real(normal[2]));
}
template< class Real >
Real Octree2< Real >::Evaluate(ConstPointer(Real) coefficients, Point3D< Real > p, const BSplineData< 2 >* fData) const
{
	Real value = Real(0);
	BSplineData< 2 > _fData;
	if (!fData) _fData.set(tree.maxDepth(), _boundaryType), fData = &_fData;
	const TreeOctNode* n = tree.nextNode();
	while (n)
	{
		Point3D< Real > c;
		Real w;
		n->centerAndWidth(c, w);
		c -= p, w *= Real(1.5);
		if (fabs(c[0]) > w || fabs(c[1]) > w || fabs(c[2]) > w)
		{
			n = tree.nextBranch(n);
			continue;
		}
		int d, off[3];
		n->depthAndOffset(d, off);
		value += (Real)
			(
				coefficients[n->nodeData.nodeIndex] *
				fData->baseFunctions[BinaryNode::CenterIndex(d, off[0])](p[0]) *
				fData->baseFunctions[BinaryNode::CenterIndex(d, off[1])](p[1]) *
				fData->baseFunctions[BinaryNode::CenterIndex(d, off[2])](p[2])
				);
		n = tree.nextNode(n);
	}
	if (_boundaryType == -1) value -= Real(0.5);
	return value;
}
template< class Real >
Pointer(Real) Octree2< Real >::Evaluate(ConstPointer(Real) coefficients, int& res, Real isoValue, int depth)
{
	int maxDepth = _boundaryType == 0 ? tree.maxDepth() - 1 : tree.maxDepth();
	if (depth <= 0 || depth > maxDepth) depth = maxDepth;
	res = 1 << depth;
	typename BSplineData< 2 >::template ValueTables< Real > vTables = _fData.template getValueTables< Real >(_fData.VALUE_FLAG);
	Pointer(Real) values = NewPointer< Real >(res * res * res);
	memset(values, 0, sizeof(Real) * res * res * res);

	for (TreeOctNode* n = tree.nextNode(); n; n = tree.nextNode(n))
	{
		if (n->depth() > (_boundaryType == 0 ? depth + 1 : depth)) continue;
		if (n->depth() < _minDepth) continue;
		int d, idx[3], start[3], end[3];
		n->depthAndOffset(d, idx);
		bool skip = false;
		for (int i = 0; i < 3; i++)
		{
			// Get the index of the functions
			idx[i] = BinaryNode::CenterIndex(d, idx[i]);
			// Figure out which samples fall into the range
			vTables.setSampleSpan(idx[i], start[i], end[i]);
			// We only care about the odd indices
			if (!(start[i] & 1)) start[i]++;
			if (!(end[i] & 1))   end[i]--;
			if (_boundaryType == 0)
			{
				// (start[i]-1)>>1 >=   res/2 
				// (  end[i]-1)<<1 <  3*res/2
				start[i] = std::max< int >(start[i], res + 1);
				end[i] = std::min< int >(end[i], 3 * res - 1);
			}
		}
		if (skip) continue;
		Real coefficient = coefficients[n->nodeData.nodeIndex];
		for (int x = start[0]; x <= end[0]; x += 2)
			for (int y = start[1]; y <= end[1]; y += 2)
				for (int z = start[2]; z <= end[2]; z += 2)
				{
					int xx = (x - 1) >> 1, yy = (y - 1) >> 1, zz = (z - 1) >> 1;
					if (_boundaryType == 0) xx -= res / 2, yy -= res / 2, zz -= res / 2;
					values[zz * res * res + yy * res + xx] +=
						coefficient *
						vTables.valueTable[idx[0] + x * vTables.functionCount] *
						vTables.valueTable[idx[1] + y * vTables.functionCount] *
						vTables.valueTable[idx[2] + z * vTables.functionCount];
				}
	}
	if (_boundaryType == -1) for (int i = 0; i < res * res * res; i++) values[i] -= Real(0.5);
	for (int i = 0; i < res * res * res; i++) values[i] -= isoValue;

	return values;
}

























////////////////
// VertexData //
////////////////
long long VertexData::CenterIndex(const TreeOctNode* node,int maxDepth)
{
	int idx[DIMENSION];
	return CenterIndex(node,maxDepth,idx);
}
long long VertexData::CenterIndex(const TreeOctNode* node,int maxDepth,int idx[DIMENSION])
{
	int d,o[3];
	node->depthAndOffset(d,o);
	for(int i=0;i<DIMENSION;i++) idx[i]=BinaryNode::CornerIndex( maxDepth+1 , d+1 , o[i]<<1 , 1 );
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
long long VertexData::CenterIndex( int depth , const int offSet[DIMENSION] , int maxDepth , int idx[DIMENSION] )
{
	for(int i=0;i<DIMENSION;i++) idx[i]=BinaryNode::CornerIndex( maxDepth+1 , depth+1 , offSet[i]<<1 , 1 );
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
long long VertexData::CornerIndex(const TreeOctNode* node,int cIndex,int maxDepth)
{
	int idx[DIMENSION];
	return CornerIndex(node,cIndex,maxDepth,idx);
}
long long VertexData::CornerIndex( const TreeOctNode* node , int cIndex , int maxDepth , int idx[DIMENSION] )
{
	int x[DIMENSION];
	Cube::FactorCornerIndex( cIndex , x[0] , x[1] , x[2] );
	int d , o[3];
	node->depthAndOffset( d , o );
	for( int i=0 ; i<DIMENSION ; i++ ) idx[i] = BinaryNode::CornerIndex( maxDepth+1 , d , o[i] , x[i] );
	return CornerIndexKey( idx );
}
long long VertexData::CornerIndex( int depth , const int offSet[DIMENSION] , int cIndex , int maxDepth , int idx[DIMENSION] )
{
	int x[DIMENSION];
	Cube::FactorCornerIndex( cIndex , x[0] , x[1] , x[2] );
	for( int i=0 ; i<DIMENSION ; i++ ) idx[i] = BinaryNode::CornerIndex( maxDepth+1 , depth , offSet[i] , x[i] );
	return CornerIndexKey( idx );
}
long long VertexData::CornerIndexKey( const int idx[DIMENSION] )
{
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
long long VertexData::FaceIndex(const TreeOctNode* node,int fIndex,int maxDepth){
	int idx[DIMENSION];
	return FaceIndex(node,fIndex,maxDepth,idx);
}
long long VertexData::FaceIndex(const TreeOctNode* node,int fIndex,int maxDepth,int idx[DIMENSION])
{
	int dir,offset;
	Cube::FactorFaceIndex(fIndex,dir,offset);
	int d,o[3];
	node->depthAndOffset(d,o);
	for(int i=0;i<DIMENSION;i++){idx[i]=BinaryNode::CornerIndex(maxDepth+1,d+1,o[i]<<1,1);}
	idx[dir]=BinaryNode::CornerIndex(maxDepth+1,d,o[dir],offset);
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
long long VertexData::EdgeIndex( const TreeOctNode* node , int eIndex , int maxDepth ){ int idx[DIMENSION] ; return EdgeIndex( node , eIndex , maxDepth , idx ); }
long long VertexData::EdgeIndex( const TreeOctNode* node , int eIndex , int maxDepth , int idx[DIMENSION] )
{
	int o , i1 , i2;
	int d , off[3];
	node->depthAndOffset( d ,off );
	Cube::FactorEdgeIndex( eIndex , o , i1 , i2 );
	for( int i=0 ; i<DIMENSION ; i++ ) idx[i] = BinaryNode::CornerIndex( maxDepth+1 , d+1 , off[i]<<1 , 1 );
	switch(o)
	{
		case 0:
			idx[1] = BinaryNode::CornerIndex( maxDepth+1 , d , off[1] , i1 );
			idx[2] = BinaryNode::CornerIndex( maxDepth+1 , d , off[2] , i2 );
			break;
		case 1:
			idx[0] = BinaryNode::CornerIndex( maxDepth+1 , d , off[0] , i1 );
			idx[2] = BinaryNode::CornerIndex( maxDepth+1 , d , off[2] , i2 );
			break;
		case 2:
			idx[0] = BinaryNode::CornerIndex( maxDepth+1 , d , off[0] , i1 );
			idx[1] = BinaryNode::CornerIndex( maxDepth+1 , d , off[1] , i2 );
			break;
	};
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
