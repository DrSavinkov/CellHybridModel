#pragma once

#include <stdlib.h>
#include <math.h>
#include <random>
#include <time.h>
#include <ctime>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <omp.h>

namespace vls {
	using namespace std;
#define DEGREE2RADIANS acos(-1.)/180.
#define Int unsigned long
	const static double PI = acos( -1. );
	const static double EPS = 1e-9;
	const static double sqr2 = pow( 2. , .5 );
	const static double sqr3 = pow( 3. , .5 );
	unsigned long long vls_longlong_pow( unsigned long long arg1 , unsigned long long arg2 ) {
		unsigned long long R = 1;
		while ( arg2 ) {
			R *= arg1;
			arg2--;
		}
		return R;
	}
	unsigned long long vls_rand( unsigned long long setmax ) {
		unsigned long long R = 0;
		int bits_used = 0;
		while ( vls_longlong_pow( 2 , bits_used ) < setmax ) {
			bits_used++;
		}
		if ( bits_used > sizeof( unsigned long long ) * 8 ) {
			printf( "Seems like setmax is too big.\n" );
			R = 0;
		}
		else {
			bool recon = true;
			while ( recon ) {
				recon = false;
				for ( int i = 1; i < bits_used; i++ ) {
					R += rand( ) % 2;
					R = R << 1;
				}
				R += rand( ) % 2;
				if ( R > setmax ) {
					R = 0;
					recon = true;
				}
			}
		}
		return R;
	}
	class PointF {
	public:
		double x , y , z;
		PointF( double a = 0. , double b = 0. , double c = 0. ) {
			x = a; y = b; z = c;
		}
		PointF( const PointF& p ) {
			x = p.x; y = p.y; z = p.z;
		}
		PointF& operator=( const PointF &p ) {
			x = p.x; y = p.y; z = p.z; return *this;
		}
		PointF operator+( const PointF &p )const {
			return PointF( x + p.x , y + p.y , z + p.z );
		}
		PointF operator-( const PointF &p )const {
			return PointF( x - p.x , y - p.y , z - p.z );
		}
		PointF operator*( const double d )const {
			return PointF( x*d , y*d , z*d );
		}
		PointF operator/( const double d )const {
			return PointF( x / d , y / d , z / d );
		}
		PointF& operator+=( const PointF &p ) {
			x += p.x; y += p.y; z += p.z; return *this;
		}
		PointF& operator-=( const PointF &p ) {
			x -= p.x; y -= p.y; z -= p.z; return *this;
		}
		PointF& operator*=( const double d ) {
			x *= d; y *= d; z *= d; return *this;
		}
		double operator|( const PointF& p ) {
			double R = x * p.x + y * p.y + z * p.z;
			return R;
		}
		PointF& operator/=( const double d ) {
			x /= d; y /= d; z /= d; return *this;
		}
		bool operator==( const PointF &p )const {
			return ( fabs( x - p.x ) < EPS ) && ( fabs( y - p.y ) < EPS ) && ( fabs( z - p.z ) < EPS );
		}
		bool operator!=( const PointF &p )const {
			return !( *this == p );
		}
		bool operator>( const PointF &p )const {
			bool r = true; if ( x < p.x ) {
				r = false;
			}
			else {
				if ( fabs( x - p.x ) < EPS ) {
					if ( y < p.y ) {
						r = false;
					}
					else {
						if ( fabs( y - p.y ) < EPS ) {
							if ( z < p.z ) {
								r = false;
							}
						}
					}
				}
			}return r;
		}
		bool operator>=( const PointF &p )const {
			return ( *this > p ) || ( *this == p );
		}
		bool operator<( const PointF &p )const {
			return !( *this >= p );
		}
		bool operator<=( const PointF &p )const {
			return !( *this > p );
		}
		char amx( )const {
			char r( 'x' ); if ( fabs( y ) > fabs( x ) )r = 'y'; if ( ( fabs( z ) > fabs( x ) ) && ( fabs( z ) > fabs( y ) ) )r = 'z'; return r;
		}
		double length( ) {
			return pow( x * x + y * y + z * z , .5 );
		}
		PointF normalize( ) {
			double L = length( );
			if ( L >= EPS ) {
				return PointF( x , y , z ) / L;
			}
			else {
				return PointF( 0 , 0 , 0 );
			}
		}
	};
	class LineF {
	public:
		PointF begin , end;
		double begin_radius , end_radius;
		LineF( PointF a = PointF( ) , PointF b = PointF( ) , double be = -1. , double en = -1. ) {
			begin = a; end = b; begin_radius = be; end_radius = en;
		}
		LineF( const LineF& l ) {
			begin = l.begin; end = l.end; begin_radius = l.begin_radius; end_radius = l.end_radius;
		}
		double length( )const {
			return pow( pow( ( begin - end ).x , 2. ) + pow( ( begin - end ).y , 2. ) + pow( ( begin - end ).z , 2. ) , .5 );
		}
		PointF mid( double d )const {
			return begin*( 1. - d ) + end*d;
		}
		vector<LineF> cut( Int n )const {
			double c( 1. / n );
			PointF pp = ( end - begin );
			vector<LineF> r;
			for ( Int i = 0; i < n; i++ )r.push_back( LineF( begin + pp*( i*c ) , begin + pp*( ( i + 1 )*c ) ) );
			return r;
		}
		PointF fzp( )const {
			return end - begin;
		}
		bool operator==( const LineF& l ) {
			return ( begin == l.begin ) && ( end == l.end );
		}
		bool operator!=( const LineF& l ) {
			return !( *this == l );
		}
		LineF& operator=( const LineF& l ) {
			begin = l.begin; end = l.end; return *this;
		}
		LineF operator+( const LineF& l )const {
			return LineF( begin , end + l.fzp( ) );
		}
		LineF operator-( const LineF& l )const {
			return LineF( begin , end - l.fzp( ) );
		}
		LineF operator*( const double d )const {
			return LineF( begin , begin + ( end - begin )*d );
		}
		LineF operator/( const double d )const {
			return LineF( begin , begin + ( end - begin ) / d );
		}
		LineF& operator+=( const LineF& l ) {
			end += l.fzp( ); return *this;
		}
		LineF& operator-=( const LineF& l ) {
			end -= l.fzp( ); return *this;
		}
		LineF& operator*=( const double d ) {
			end = begin + ( end - begin )*d; return *this;
		}
		LineF& operator/=( const double d ) {
			end = begin + ( end - begin ) / d; return *this;
		}
		double operator*( const LineF& l )const {
			return fzp( ).x*l.fzp( ).x + fzp( ).y*l.fzp( ).y + fzp( ).z*l.fzp( ).z;
		}
		LineF orth( )const {
			PointF ao; if ( length( ) < EPS ) {
				ao = PointF( 1. );
			}
			else {
				PointF wv = fzp( );/* (a, b, c), a != 0 --> (e, f, k), f = 1 --> (e, 1, k) --> ae+b+ck = 0 --> k = 1 --> (e, 1, 1) --> ae + b + c = 0 --> e = (-b-c)/a */if ( wv.amx( ) == 'x' )ao = PointF( ( -wv.y - wv.z ) / wv.x , 1 , 1 ); if ( wv.amx( ) == 'y' )ao = PointF( 1 , ( -wv.x - wv.z ) / wv.y , 1 ); if ( wv.amx( ) == 'z' )ao = PointF( 1 , 1 , ( -wv.x - wv.y ) / wv.z ); ao /= LineF( PointF( ) , ao ).length( );
			}return LineF( begin , begin + ao );
		}
		LineF rot( LineF rotor = LineF( PointF( ) , PointF( 1. ) ) , double degree = 0. )const {
			LineF r; if ( rotor.length( ) < EPS ) {
				rotor = LineF( PointF( ) , PointF( 1. ) ); degree = 0.;
			}double radians = DEGREE2RADIANS*degree; if ( length( ) < EPS ) {
				r = *this;
			}
			else {
				PointF prot = rotor.fzp( ) / rotor.length( ); PointF pvec = this->fzp( ) / this->length( ); double i = pvec.x*( cos( radians ) + ( 1. - cos( radians ) )*prot.x*prot.x ) + pvec.y*( ( 1. - cos( radians ) )*prot.x*prot.y - sin( radians )*prot.z ) + pvec.z*( ( 1. - cos( radians ) )*prot.x*prot.z + sin( radians )*prot.y ); double j = pvec.x*( ( 1. - cos( radians ) )*prot.y*prot.x + sin( radians )*prot.z ) + pvec.y*( cos( radians ) + ( 1. - cos( radians ) )*prot.y*prot.y ) + pvec.z*( ( 1. - cos( radians ) )*prot.y*prot.z - sin( radians )*prot.x ); double k = pvec.x*( ( 1. - cos( radians ) )*prot.z*prot.x - sin( radians )*prot.y ) + pvec.y*( ( 1. - cos( radians ) )*prot.z*prot.y + sin( radians )*prot.x ) + pvec.z*( cos( radians ) + ( 1. - cos( radians ) )*prot.z*prot.z ); r = LineF( begin , begin + PointF( i , j , k )*length( ) );
			}return r;
		}
		LineF tsl( ) {
			return *this / this->length( );
		}
	};
	class TreeF {
	public:
		double maxLength;
		vector<PointF> points;
		map<Int/*start*/ , map<Int/*finish*/ , double/*length*/> > connections;
		void insert_line( const LineF& l , double length = -1. ) {
			map<PointF , Int> fastLinks;
			long long begin_id = -1 , end_id = -1;
			for ( Int i = 0; i < points.size( ); i++ ) {
				if ( points [ i ] == l.begin ) {
					begin_id = i;
				}
				if ( points [ i ] == l.end ) {
					end_id = i;
				}
			}
			if ( length < -.5 )length = l.length( );
			if ( length > maxLength ) {
				double divTo = length / maxLength + 1.;
				vector<LineF> subsectors = l.cut( Int( divTo ) );
				for ( Int i = 0; i < subsectors.size( ); i++ ) {
					insert_line( subsectors [ i ] , length / Int( divTo ) );
				}
			}
			else {
				if ( begin_id == -1 ) {
					begin_id = points.size( );
					points.push_back( l.begin );
				}
				if ( end_id == -1 ) {
					end_id = points.size( );
					points.push_back( l.end );
				}
				connections [ unsigned long( begin_id ) ] [ unsigned long( end_id ) ] = length;
				connections [ unsigned long( end_id ) ] [ unsigned long( begin_id ) ] = length;
			}
		}
		void addTree( TreeF &t ) {

		}
		TreeF( ) {
			maxLength = 1.;
		}
		TreeF& operator=( const TreeF& t ) {
			maxLength = t.maxLength; points = t.points; connections = t.connections; return *this;
		}
		void throwToScriptJustNodes( string fileName , double coefficiento ) {
			map<Int , map<Int , double> > justNodes;
			for ( Int i = 0; i < points.size( ); i++ ) {
				if ( connections [ i ].size( ) != 2 ) {
					vector<Int> foundPoints;
					for ( auto j = connections [ i ].begin( ); j != connections [ i ].end( ); j++ ) {
						auto prevpoint = i , currpoint = j->first;
						while ( connections [ currpoint ].size( ) == 2 ) {
							auto buf = prevpoint;
							prevpoint = currpoint;
							auto link = connections [ currpoint ].begin( );
							if ( link->first == buf ) {
								link++;
							}
							currpoint = link->first;
						}
						foundPoints.push_back( currpoint );
					}
					for ( Int j = 0; j < foundPoints.size( ); j++ ) {
						justNodes [ i ] [ foundPoints [ j ] ] = LineF( points [ i ] , points [ foundPoints [ j ] ] ).length( );
					}
				}
			}
			string ssx;
			FILE* foutM;
			double swapVecMassive [ 6 ];
			ssx = fileName + string( ".obj" );
			foutM = fopen( ssx.c_str( ) , "w" );
			for ( map<Int , map<Int , double> >::iterator i = justNodes.begin( ); i != justNodes.end( ); i++ ) {
				for ( map<Int , double>::iterator j = i->second.begin( ); j != i->second.end( ); j++ ) {
					swapVecMassive [ 0 ] = points [ i->first ].x;
					swapVecMassive [ 1 ] = points [ i->first ].y;
					swapVecMassive [ 2 ] = points [ i->first ].z;
					swapVecMassive [ 3 ] = points [ j->first ].x;
					swapVecMassive [ 4 ] = points [ j->first ].y;
					swapVecMassive [ 5 ] = points [ j->first ].z;
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 0 ] + coefficiento / 2.0 ,
							 swapVecMassive [ 1 ] ,
							 swapVecMassive [ 2 ]
					);
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 3 ] + coefficiento / 2.0 ,
							 swapVecMassive [ 4 ] ,
							 swapVecMassive [ 5 ]
					);
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 3 ] - coefficiento / 2.0 ,
							 swapVecMassive [ 4 ] ,
							 swapVecMassive [ 5 ]
					);
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 0 ] - coefficiento / 2.0 ,
							 swapVecMassive [ 1 ] ,
							 swapVecMassive [ 2 ]
					);
					///==============================================
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 0 ] ,
							 swapVecMassive [ 1 ] + coefficiento / 2.0 ,
							 swapVecMassive [ 2 ]
					);
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 3 ] ,
							 swapVecMassive [ 4 ] + coefficiento / 2.0 ,
							 swapVecMassive [ 5 ]
					);
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 3 ] ,
							 swapVecMassive [ 4 ] - coefficiento / 2.0 ,
							 swapVecMassive [ 5 ]
					);
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 0 ] ,
							 swapVecMassive [ 1 ] - coefficiento / 2.0 ,
							 swapVecMassive [ 2 ]
					);
					///==============================================
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 0 ] ,
							 swapVecMassive [ 1 ] ,
							 swapVecMassive [ 2 ] + coefficiento / 2.0
					);
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 3 ] ,
							 swapVecMassive [ 4 ] ,
							 swapVecMassive [ 5 ] + coefficiento / 2.0
					);
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 3 ] ,
							 swapVecMassive [ 4 ] ,
							 swapVecMassive [ 5 ] - coefficiento / 2.0
					);
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 0 ] ,
							 swapVecMassive [ 1 ] ,
							 swapVecMassive [ 2 ] - coefficiento / 2.0
					);
				}
			}
			long i = 0;
			for ( map<Int , map<Int , double> >::iterator ii = justNodes.begin( ); ii != justNodes.end( ); ii++ ) {
				for ( map<Int , double>::iterator j = ii->second.begin( ); j != ii->second.end( ); j++ ) {
					fprintf( foutM , "f %ld %ld %ld\nf %ld %ld %ld\nf %ld %ld %ld\nf %ld %ld %ld\nf %ld %ld %ld\nf %ld %ld %ld\n" ,
							 i * 6 + 1 , i * 6 + 2 , i * 6 + 3 ,
							 i * 6 + 3 , i * 6 + 4 , i * 6 + 1 ,
							 i * 6 + 5 , i * 6 + 6 , i * 6 + 7 ,
							 i * 6 + 7 , i * 6 + 8 , i * 6 + 5 ,
							 i * 6 + 9 , i * 6 + 10 , i * 6 + 11 ,
							 i * 6 + 11 , i * 6 + 12 , i * 6 + 9
					);
					i += 2;
				}
			}
			fclose( foutM );
		}
		void throwToScript( string fileName , double coefficiento ) {
			string ssx;
			FILE* foutM;
			double swapVecMassive [ 6 ];
			ssx = fileName + string( ".obj" );
			foutM = fopen( ssx.c_str( ) , "w" );
			for ( map<Int , map<Int , double> >::iterator i = connections.begin( ); i != connections.end( ); i++ ) {
				for ( map<Int , double>::iterator j = i->second.begin( ); j != i->second.end( ); j++ ) {
					swapVecMassive [ 0 ] = points [ i->first ].x;
					swapVecMassive [ 1 ] = points [ i->first ].y;
					swapVecMassive [ 2 ] = points [ i->first ].z;
					swapVecMassive [ 3 ] = points [ j->first ].x;
					swapVecMassive [ 4 ] = points [ j->first ].y;
					swapVecMassive [ 5 ] = points [ j->first ].z;
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 0 ] + coefficiento / 2.0 ,
							 swapVecMassive [ 1 ] ,
							 swapVecMassive [ 2 ]
					);
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 3 ] + coefficiento / 2.0 ,
							 swapVecMassive [ 4 ] ,
							 swapVecMassive [ 5 ]
					);
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 3 ] - coefficiento / 2.0 ,
							 swapVecMassive [ 4 ] ,
							 swapVecMassive [ 5 ]
					);
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 0 ] - coefficiento / 2.0 ,
							 swapVecMassive [ 1 ] ,
							 swapVecMassive [ 2 ]
					);
					///==============================================
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 0 ] ,
							 swapVecMassive [ 1 ] + coefficiento / 2.0 ,
							 swapVecMassive [ 2 ]
					);
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 3 ] ,
							 swapVecMassive [ 4 ] + coefficiento / 2.0 ,
							 swapVecMassive [ 5 ]
					);
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 3 ] ,
							 swapVecMassive [ 4 ] - coefficiento / 2.0 ,
							 swapVecMassive [ 5 ]
					);
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 0 ] ,
							 swapVecMassive [ 1 ] - coefficiento / 2.0 ,
							 swapVecMassive [ 2 ]
					);
					///==============================================
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 0 ] ,
							 swapVecMassive [ 1 ] ,
							 swapVecMassive [ 2 ] + coefficiento / 2.0
					);
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 3 ] ,
							 swapVecMassive [ 4 ] ,
							 swapVecMassive [ 5 ] + coefficiento / 2.0
					);
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 3 ] ,
							 swapVecMassive [ 4 ] ,
							 swapVecMassive [ 5 ] - coefficiento / 2.0
					);
					fprintf( foutM , "v %lf %lf %lf\n" ,
							 swapVecMassive [ 0 ] ,
							 swapVecMassive [ 1 ] ,
							 swapVecMassive [ 2 ] - coefficiento / 2.0
					);
				}
			}
			long i = 0;
			for ( map<Int , map<Int , double> >::iterator ii = connections.begin( ); ii != connections.end( ); ii++ ) {
				for ( map<Int , double>::iterator j = ii->second.begin( ); j != ii->second.end( ); j++ ) {
					fprintf( foutM , "f %ld %ld %ld\nf %ld %ld %ld\nf %ld %ld %ld\nf %ld %ld %ld\nf %ld %ld %ld\nf %ld %ld %ld\n" ,
							 i * 6 + 1 , i * 6 + 2 , i * 6 + 3 ,
							 i * 6 + 3 , i * 6 + 4 , i * 6 + 1 ,
							 i * 6 + 5 , i * 6 + 6 , i * 6 + 7 ,
							 i * 6 + 7 , i * 6 + 8 , i * 6 + 5 ,
							 i * 6 + 9 , i * 6 + 10 , i * 6 + 11 ,
							 i * 6 + 11 , i * 6 + 12 , i * 6 + 9
					);
					i += 2;
				}
			}
			fclose( foutM );
		}
		void throwToPOBJ( string filename ) {
			FILE* f = fopen( ( filename + ".pobj" ).c_str( ) , "w" );
			if ( f ) {
				//fprintf(f, "POINTS %ld\n", points.size());
				Int lines = 0;
				for ( auto i = connections.begin( ); i != connections.end( ); i++ ) {
					lines += unsigned long( i->second.size( ) );
				}
				//fprintf(f, "LINES %ld\n", lines);
				for ( Int i = 0; i < points.size( ); i++ ) {
					fprintf( f , "P %lf %lf %lf\n" , points [ i ].x , points [ i ].y , points [ i ].z );
				}
				for ( auto i = connections.begin( ); i != connections.end( ); i++ ) {
					for ( auto j = i->second.begin( ); j != i->second.end( ); j++ ) {
						fprintf( f , "L %ld %ld\n" , i->first + 1 , j->first + 1 );
					}
				}
				fclose( f );
			}
		}
		map<unsigned long , map<unsigned long , double> > getJustNodes( ) {
			map<unsigned long , map<unsigned long , double> > justNodes;
			for ( unsigned long i = 0; i < points.size( ); i++ ) {
				if ( connections [ i ].size( ) > 2 ) {
					vector<unsigned long> foundPoints;
					vector<double> foundlength;
					for ( auto j = connections [ i ].begin( ); j != connections [ i ].end( ); j++ ) {
						auto prevpoint = i , currpoint = j->first;
						double reallength = LineF( points [ prevpoint ] , points [ currpoint ] ).length( );
						while ( connections [ currpoint ].size( ) == 2 ) {
							auto buf = prevpoint;
							prevpoint = currpoint;
							auto link = connections [ currpoint ].begin( );
							if ( link->first == buf ) {
								link++;
							}
							currpoint = link->first;
							reallength += LineF( points [ prevpoint ] , points [ currpoint ] ).length( );
						}
						if ( connections [ currpoint ].size( ) > 2 ) {
							foundPoints.push_back( currpoint );
							foundlength.push_back( reallength );
						}
					}
					for ( unsigned long j = 0; j < foundPoints.size( ); j++ ) {
						justNodes [ i ] [ foundPoints [ j ] ] = foundlength [ j ];// LineF(points[i], points[foundPoints[j]]).length();
					}
				}
			}
			return justNodes;
		}
		void getFromPOBJ( string file ) {
			FILE* fin = fopen( file.c_str( ) , "r" );
			char read [ 4096 ];
			string sread;
			if ( fin ) {
				while ( fscanf( fin , "%s\0" , read ) > 0 ) {
					sread = string( read );
					if ( sread == "P" || sread == "p" ) {
						double x , y , z;
						fscanf( fin , "%s\0" , read );
						sscanf( read , "%lf" , &x );
						fscanf( fin , "%s\0" , read );
						sscanf( read , "%lf" , &y );
						fscanf( fin , "%s\0" , read );
						sscanf( read , "%lf" , &z );
						points.push_back( PointF( x , y , z ) );
					}
					if ( sread == "L" || sread == "l" ) {
						unsigned long gin , gout;
						fscanf( fin , "%s\0" , read );
						sscanf( read , "%ld" , &gin );
						fscanf( fin , "%s\0" , read );
						sscanf( read , "%ld" , &gout );
						gin--;
						gout--;
						connections [ gin ] [ gout ] = ( points [ gin ] - points [ gout ] ).length( );
					}
				}
				fclose( fin );
			}
		}
	};
	double rnd( bool positive = false );
	void write( const PointF& p , bool el = true ) {
		printf( "PF (%lf,%lf,%lf)" , p.x , p.y , p.z ); if ( el )printf( "\n" );
	}
	void write( const LineF& l , bool el = true ) {
		printf( "LF (%lf,%lf,%lf)-->(%lf,%lf,%lf)" , l.begin.x , l.begin.y , l.begin.z , l.end.x , l.end.y , l.end.z ); if ( el )printf( "\n" );
	}
	void write( const TreeF& t , bool el = true ) {
		for ( Int i = 0; i < t.points.size( ); i++ ) {
			write( t.points [ i ] , el );
		}
	}
	void write( string s = "" , bool el = true ) {
		if ( el )printf( "%s\n" , s.c_str( ) ); if ( !el )printf( "%s" , s.c_str( ) );
	}
	void write( double v , bool el = true ) {
		if ( !el )printf( "%lf" , v ); if ( el )printf( "%lf\n" , v );
	}
	void write( Int v , bool el = true ) {
		if ( !el )printf( "%ld" , v ); if ( el )printf( "%ld\n" , v );
	}
	void write( long v , bool el = true ) {
		if ( !el )printf( "%ld" , v ); if ( el )printf( "%ld\n" , v );
	}
	void write( int v , bool el = true ) {
		if ( !el )printf( "%ld" , v ); if ( el )printf( "%ld\n" , v );
	}
	void write( unsigned int v , bool el = true ) {
		if ( !el )printf( "%i" , v ); if ( el )printf( "%i\n" , v );
	}
	string sprint( double v ) {
		char r [ 4096 ]; sprintf( r , "%lf\0" , v ); return string( r );
	}
	string sprint( int v ) {
		char r [ 4096 ]; sprintf( r , "%d\0" , v ); return string( r );
	}
	string sprint( long v ) {
		char r [ 4096 ]; sprintf( r , "%ld\0" , v ); return string( r );
	}
	string sprint( Int v ) {
		char r [ 4096 ]; sprintf( r , "%ld\0" , v ); return string( r );
	}
	void wait( Int msec ) {
		double wtime = msec*.001; double startTime = omp_get_wtime( ); while ( omp_get_wtime( ) - startTime < wtime ) {
		}
	}
	double rnd( bool positive ) {
		//double r( .5*( 1. - 2.*( rand( ) % 2 ) )*( sin( rand( ) ) + cos( rand( ) ) ) ); if (r < 0.&&positive)r = -r; return r;
		double r = ( double( rand( ) ) / double( RAND_MAX ) ) * ( 1 - 2 * ( rand( ) % 2 ) ); if ( positive && ( r < 0. ) ) {
			r = -r;
		}
		return r;
	}
	double gravityIteration(
		TreeF& t ,
		map<long , map<long , map<long , vector<Int> > > >& zones ,
		map<long , map<long , map<long , vector<Int> > > >& zonesWork ,
		double wantNCDistance = 1. ,
		bool divSum = true ,
		map<long , map<long , map<long , PointF> > > &randVecFlow = map<long , map<long , map<long , PointF> > >( ) ,
		double dt = 1.0 ,
		double nomi = 0.0 ,
		double zoneSizeGet = 10. ,
		bool *zoneCollisionsGood = NULL ,
		set<PointF> &staticPoints = set<PointF>( ) ,
		vector<PointF> &speedSave = vector<PointF>( ) ,
		double speedSaving = 0.0 ,
		int *collized = NULL ,
		double farFromBoard = 0. ,
		bool ignoreFarPoints = false ,
		double *ad = NULL
	) {
		int adcount = 0;
		if ( collized ) {
			*collized = 0;
		}
		const double coeffUse = .5;
		double WNCD = 0;
		WNCD = wantNCDistance;
		const double zoneSize = zoneSizeGet;
		for ( auto i = zones.begin( ); i != zones.end( ); i++ ) {
			for ( auto j = i->second.begin( ); j != i->second.end( ); j++ ) {
				for ( auto k = j->second.begin( ); k != j->second.end( ); k++ ) {
					k->second.clear( );
				}
			}
		}
		for ( long i = 0; i < t.points.size( ); i++ ) {
			long zx = long( t.points [ i ].x / zoneSize + int( t.points [ i ].x > ( zoneSize*.5 ) )*.5 - int( t.points [ i ].x < ( -zoneSize*.5 ) )*.5 );
			long zy = long( t.points [ i ].y / zoneSize + int( t.points [ i ].y > ( zoneSize*.5 ) )*.5 - int( t.points [ i ].y < ( -zoneSize*.5 ) )*.5 );
			long zz = long( t.points [ i ].z / zoneSize + int( t.points [ i ].z > ( zoneSize*.5 ) )*.5 - int( t.points [ i ].z < ( -zoneSize*.5 ) )*.5 );
			zones [ zx ] [ zy ] [ zz ].push_back( i );
		}
		for ( auto i = zonesWork.begin( ); i != zonesWork.end( ); i++ ) {
			for ( auto j = i->second.begin( ); j != i->second.end( ); j++ ) {
				for ( auto k = j->second.begin( ); k != j->second.end( ); k++ ) {
					k->second.clear( );
				}
			}
		}
		map<long , map<long , map<long , vector<Int> > > >::iterator iter_1;
		map<long , map<long , vector<Int> > >::iterator iter_2;
		map<long , vector<Int> >::iterator iter_3;
		for ( auto i = zones.begin( ); i != zones.end( ); i++ ) {
			for ( auto j = i->second.begin( ); j != i->second.end( ); j++ ) {
				for ( auto k = j->second.begin( ); k != j->second.end( ); k++ ) {
					for ( int dx = -1; dx < 2; dx++ ) {
						for ( int dy = -1; dy < 2; dy++ ) {
							for ( int dz = -1; dz < 2; dz++ ) {
								iter_1 = zones.find( i->first + dx );
								if ( iter_1 != zones.end( ) ) {
									iter_2 = iter_1->second.find( j->first + dy );
									if ( iter_2 != iter_1->second.end( ) ) {
										iter_3 = iter_2->second.find( k->first + dz );
										if ( iter_3 != iter_2->second.end( ) ) {
											for ( Int f = 0; f < zones [ i->first + dx ] [ j->first + dy ] [ k->first + dz ].size( ); f++ ) {
												zonesWork [ i->first ] [ j->first ] [ k->first ].push_back( zones [ i->first + dx ] [ j->first + dy ] [ k->first + dz ] [ f ] );
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		double diff = 0.;
		vector<PointF> dim [ 2 ] = { t.points,t.points };
		double mdc = 1. / t.points.size( );
		vector<PointF> deltas;
		deltas.resize( dim [ 0 ].size( ) );
		if ( speedSave.size( ) < dim [ 0 ].size( ) ) {
			speedSave.resize( dim [ 0 ].size( ) );
#pragma omp parallel for
			for ( long i = 0; i < dim [ 0 ].size( ); i++ ) {
				speedSave [ i ] = PointF( );
			}
		}
#pragma omp parallel for
		for ( long i = 0; i < dim [ 0 ].size( ); i++ ) {
			deltas [ i ] = PointF( );
		}
		if ( zoneCollisionsGood != NULL ) {
			*zoneCollisionsGood = true;
		}
#pragma omp parallel for
		for ( long i = 0; i < dim [ 0 ].size( ); i++ ) {
			Int idI = i;
			PointF delta( 0. , 0. , 0. );
			PointF delta2 = delta;
			double counter = 0;
			long zx = long( dim [ 0 ] [ i ].x / zoneSize + int( dim [ 0 ] [ i ].x > ( zoneSize*.5 ) )*.5 - int( dim [ 0 ] [ i ].x < ( -zoneSize*.5 ) )*.5 );
			long zy = long( dim [ 0 ] [ i ].y / zoneSize + int( dim [ 0 ] [ i ].y > ( zoneSize*.5 ) )*.5 - int( dim [ 0 ] [ i ].y < ( -zoneSize*.5 ) )*.5 );
			long zz = long( dim [ 0 ] [ i ].z / zoneSize + int( dim [ 0 ] [ i ].z > ( zoneSize*.5 ) )*.5 - int( dim [ 0 ] [ i ].z < ( -zoneSize*.5 ) )*.5 );
			if ( randVecFlow.find( zx ) != randVecFlow.end( ) ) {
				if ( randVecFlow [ zx ].find( zy ) != randVecFlow [ zx ].end( ) ) {
					if ( randVecFlow [ zx ] [ zy ].find( zz ) != randVecFlow [ zx ] [ zy ].end( ) ) {
						PointF lq = PointF( zx , zy , zz ) * zoneSize;
						double fds = farFromBoard - LineF( lq , dim [ 0 ] [ i ] ).length( );
						delta2 += randVecFlow [ zx ] [ zy ] [ zz ] * ( fds * int( fds > 0. ) + 1. );
						if ( zoneCollisionsGood != NULL && ( farFromBoard > 0. ) ) {
							*zoneCollisionsGood = false;
						}
						if ( collized != NULL && ( farFromBoard > 0. ) ) {
#pragma omp atomic
							*collized += 1;
						}
					}
				}
			}
			vector<Int> lci;
			iter_1 = zonesWork.find( zx );
			if ( iter_1 != zonesWork.end( ) ) {
				iter_2 = iter_1->second.find( zy );
				if ( iter_2 != iter_1->second.end( ) ) {
					iter_3 = iter_2->second.find( zz );
					if ( iter_3 != iter_2->second.end( ) ) {
						lci = zonesWork [ zx ] [ zy ] [ zz ];
					}
				}
			}
			double maxDeltaVec = 0.;
			for ( Int j = 0; j < lci.size( ); j++ ) {
				if ( i != lci [ j ] ) {
					Int idJ = lci [ j ];
					double realDistance = LineF( dim [ 0 ] [ idI ] , dim [ 0 ] [ idJ ] ).length( );
					bool moveClose = ( t.connections [ idI ].find( idJ ) != t.connections [ idI ].end( ) );
					if ( !moveClose ) {
						if ( realDistance < WNCD ) {
							// растолкнуть
							bool tooCloseInGraph = false;
							for ( auto tcig = t.connections [ idI ].begin( ); tcig != t.connections [ idI ].end( ); tcig++ ) {
								if ( t.connections [ tcig->first ].find( idJ ) != t.connections [ tcig->first ].end( ) ) {
									tooCloseInGraph = true;
								}
							}
							if ( tooCloseInGraph || ( !ignoreFarPoints ) ) {
								if ( ad != NULL ) {
#pragma omp atomic
									*ad += WNCD - realDistance;
									adcount++;
								}
								PointF DL = LineF( dim [ 0 ] [ idJ ] , dim [ 0 ] [ idI ] ).tsl( ).fzp( )*coeffUse*( WNCD - realDistance );
								double DL_d = DL.length( );
								if ( maxDeltaVec < DL_d ) {
									maxDeltaVec = DL_d;
								}
								delta += DL;
								counter++;
							}
						}
					}
				}
			}
			PointF pdelta( 0 , 0 , 0 );
			if ( t.connections.find( idI ) != t.connections.end( ) ) {
				for ( map<Int , double>::iterator j = t.connections [ idI ].begin( ); j != t.connections [ idI ].end( ); j++ ) {
					Int idJ = j->first;
					double realDistance = LineF( dim [ 0 ] [ idI ] , dim [ 0 ] [ idJ ] ).length( );
					bool moveClose = ( t.connections [ idI ].find( idJ ) != t.connections [ idI ].end( ) );
					if ( moveClose ) {
						double trgL = t.connections [ idI ] [ idJ ];
						// отправить на правильные места
						PointF DL = LineF( dim [ 0 ] [ idI ] , dim [ 0 ] [ idJ ] ).tsl( ).fzp( )*coeffUse*( realDistance - trgL );
						double DL_d = DL.length( );
						if ( maxDeltaVec < DL_d ) {
							maxDeltaVec = DL_d;
						}
						delta += DL;
						counter++;
					}
				}
			}
			if ( counter > 0 )delta /= counter;
			PointF demiP = PointF( rnd( ) , rnd( ) , rnd( ) );
			demiP = LineF( PointF( ) , demiP ).tsl( ).fzp( ) * nomi;
			demiP = LineF( PointF( ) , dim [ 0 ] [ i ] ).tsl( ).fzp( ) * nomi;
			delta += demiP;
			delta = delta * dt;
			if ( staticPoints.find( dim [ 0 ] [ i ] ) != staticPoints.end( ) ) {
				delta2 = delta = PointF( );
				dim [ 1 ] [ i ] = dim [ 0 ] [ i ];
			}
			else {
				dim [ 1 ] [ i ] = dim [ 0 ] [ i ] + delta + ( speedSave [ i ] * dt ) + delta2 * dt;
				speedSave [ i ] += ( delta / dt ) + ( delta2 / dt );
				speedSave [ i ] *= speedSaving;
			}
		}
		if ( adcount ) {
			*ad /= adcount;
		}
		for ( long i = 0; i < dim [ 1 ].size( ); i++ ) {
			Int idI = i;
			if ( t.connections.find( idI ) != t.connections.end( ) ) {
				for ( auto j = t.connections [ idI ].begin( ); j != t.connections [ idI ].end( ); j++ ) {
					Int idJ = j->first;
					double realDistance = LineF( dim [ 1 ] [ idI ] , dim [ 1 ] [ idJ ] ).length( );
					bool moveClose = ( t.connections [ idI ].find( idJ ) != t.connections [ idI ].end( ) );
					if ( moveClose ) {
						double trgL = t.connections [ idI ] [ idJ ];
						diff += fabs( realDistance - trgL );
					}
				}
			}
		}
		t.points = dim [ 1 ];
		if ( divSum ) {
			diff /= dim [ 1 ].size( );
		}
		return diff;
	}
}