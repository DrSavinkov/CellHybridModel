// CellHybridModel.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "VLS.h"

using namespace vls;

class lmtr {
public:
	long a , b , c;
};
vector<PointF> lmpoints;
vector<vector<lmtr>> lmtries;
vector<PointF> lmColors;
class lmModel {
public:
	vector<PointF> lmpoints;
	vector<lmtr> lmtries;
	void loadModel( string s ) {
		ifstream fin( s , ifstream::in );
		if ( fin.is_open( ) ) {
			lmtries.clear( );
			//cout << "Fin \"" << s << "\" is open." << endl;
			while ( !fin.eof( ) ) {
				fin >> s;
				if ( s == "v" ) {
					double x , y , z;
					fin >> x >> y >> z;
					lmpoints.push_back( PointF( x , y , z ) );
				}
				if ( s == "f" ) {
					long x , y , z;
					fin >> x >> y >> z;
					x--;
					y--;
					z--;
					lmtr b;
					b.a = x;
					b.b = y;
					b.c = z;
					lmtries.push_back( b );
				}
			}
			fin.close( );
		}
	}
};
void loadModelForRezult( string s , PointF POS , double zoom , PointF color ) {
	long pid = long( lmpoints.size( ) );
	ifstream fin( s , ifstream::in );
	if ( fin.is_open( ) ) {
		vector<lmtr> tries;
		//cout << "Fin \"" << s << "\" is open." << endl;
		while ( !fin.eof( ) ) {
			fin >> s;
			if ( s == "v" ) {
				double x , y , z;
				fin >> x >> y >> z;
				lmpoints.push_back( PointF( x*zoom + POS.x , y*zoom + POS.y , z*zoom + POS.z ) );
			}
			if ( s == "f" ) {
				long x , y , z;
				fin >> x >> y >> z;
				x--;
				y--;
				z--;
				lmtr b;
				b.a = x + pid;
				b.b = y + pid;
				b.c = z + pid;
				tries.push_back( b );
			}
		}
		if ( tries.size( ) ) {
			lmColors.push_back( color );
			lmtries.push_back( tries );
		}
		fin.close( );
	}
}
void pushLmModel( lmModel md , PointF POS , double zoom , PointF color ) {
	long pid = long( lmpoints.size( ) );
	vector<lmtr> tries = md.lmtries;
	for ( long i = 0; i < md.lmpoints.size( ); i++ ) {
		lmpoints.push_back( md.lmpoints [ i ] * zoom + POS );
	}
	for ( long i = 0; i < tries.size( ); i++ ) {
		tries [ i ].a += pid;
		tries [ i ].b += pid;
		tries [ i ].c += pid;
	}
	if ( tries.size( ) ) {
		lmColors.push_back( color );
		lmtries.push_back( tries );
	}
}
void throwAllModels2Out( string s ) {
	ofstream fout( s , ofstream::out );
	if ( fout.is_open( ) ) {
		for ( long i = 0; i < lmpoints.size( ); i++ ) {
			fout << "v " << lmpoints [ i ].x << " " << lmpoints [ i ].y << " " << lmpoints [ i ].z << endl;
		}
		for ( long i = 0; i < lmtries.size( ); i++ ) {
			fout << "#_NEW_MODEL " << lmColors [ i ].x << " " << lmColors [ i ].y << " " << lmColors [ i ].z << endl;
			for ( long j = 0; j < lmtries [ i ].size( ); j++ ) {
				fout << "f " << lmtries [ i ] [ j ].a + 1 << " " << lmtries [ i ] [ j ].b + 1 << " " << lmtries [ i ] [ j ].c + 1 << endl;
			}
		}
		fout.close( );
	}
}

void wait( double seconds ) {
	double ct = omp_get_wtime( );
	while ( omp_get_wtime( ) - ct < seconds );
}
enum FLUID_ZONE_TYPE {
	NORMAL , BOUND , FRC , VESSEL , B_CELL
};
enum CELL_ID {
	ID_CD4p , ID_CD4pi , ID_CD8p , ID_OHTER
};
vector<double> GTC; // GLOBAL TIME COUNTER
void gtc_start( ) {
	GTC.push_back( omp_get_wtime( ) );
}
double gtc_stop( ) {
	double R = 0;
	if ( GTC.size( ) > 0 ) {
		R = omp_get_wtime( ) - GTC [ GTC.size( ) - 1 ];
	}
	return R;
}
PointF noize( ) {
	return PointF( rnd( ) , rnd( ) , rnd( ) ).normalize( ) * EPS;
}
const static long FLUID_ZONE_N = 1;
#define GROWTH_FACTOR 0
#define INFLAMMATORY_FACTOR 1
#define VIRUS_FACTOR 2
const static double ft = 1.; // 1 hour
const static double dt = 0.1; // discretization 0.05*0.05 / 10. = 0.00025
const static double MAX_CELL_DIAMETER = 6.;
const static double VEC_MOD = 0.75; // доля старого вектора, оставляемая на итерации
const static double R_MOD = 0.1; // доля вектора, определяемая случайным значением
const static double WALL__TOUCH_FORCE_COEFFICIENT = 0.6 * ( 1. - VEC_MOD - R_MOD * ( 1. - VEC_MOD ) );
const static double CELL__TOUCH_FORCE_COEFFICIENT = 0.2 * ( 1. - VEC_MOD - R_MOD * ( 1. - VEC_MOD ) );
const static double FLUID_TOUCH_FORCE_COEFFICIENT = 0.3 * ( 1. - VEC_MOD - R_MOD * ( 1. - VEC_MOD ) );
const static long ITERS = long( ( ft / dt ) + 0.5 );
double degradation [ FLUID_ZONE_N ] = { 0. };
double RND( ) {
	return ( double( rand( ) ) ) / RAND_MAX;
}
class PointL {
public:
	long x , y , z;
	PointL( PointF p = PointF( ) ) {
		x = long( p.x + ( int( p.x > 0 ) - int( p.x < 0 ) ) * 0.5 );
		y = long( p.y + ( int( p.y > 0 ) - int( p.y < 0 ) ) * 0.5 );
		z = long( p.z + ( int( p.z > 0 ) - int( p.z < 0 ) ) * 0.5 );
	}
	PointL( double X , double Y , double Z ) {
		x = long( X + ( int( X > 0 ) - int( X < 0 ) ) * 0.5 );
		y = long( Y + ( int( Y > 0 ) - int( Y < 0 ) ) * 0.5 );
		z = long( Z + ( int( Z > 0 ) - int( Z < 0 ) ) * 0.5 );
	}
	PointF convert( ) {
		return PointF( x , y , z );
	}
};
class CELL {
public:
	CELL_ID ID;
	PointF POS;
	PointF VEC;
	double LIFE , LIFE_CONST;
	double DIV_INTERVAL , DIV_INTERVAL_CONST;
	double MITOS , MITOS_CONST;
	double SPEED;
	double RADIUS;
	bool mitos( ) {
		return ( MITOS > 0. );
	}
	bool live( ) {
		return ( LIFE > 0 );
	}
	bool can_div( ) {
		return ( DIV_INTERVAL <= 0. ) && ( MITOS <= 0 ) && live( );
	}
	bool dead( ) {
		return ( LIFE <= 0 );
	}
	void stat_time_update( double dt , bool& create_new_cell ) {
		create_new_cell = ( ( MITOS > 0 ) && ( MITOS <= dt ) );
		LIFE -= dt;
		DIV_INTERVAL -= dt;
		MITOS -= dt;
		if ( LIFE < 0 ) {
			LIFE = 0;
		}
		if ( DIV_INTERVAL < 0 ) {
			DIV_INTERVAL = 0;
		}
		if ( MITOS < 0 ) {
			MITOS = 0;
		}
	}
	bool can_move( ) {
		return MITOS == 0;
	}
};
class CD4p : public CELL {
public:
	CD4p( PointF P = PointF( ) , CELL_ID id = ID_CD4p ) {
		VEC = PointF( );
		POS = P;
		DIV_INTERVAL = 0;
		MITOS = 0;
		LIFE = 14400 * int( id == ID_CD4p ) + 2160 * int( id == ID_CD4pi );
		DIV_INTERVAL_CONST = 300;
		MITOS_CONST = 60;
		LIFE_CONST = 14400 * int( id == ID_CD4p ) + 2160 * int( id == ID_CD4pi );
		SPEED = 10;
		RADIUS = 3.;
		ID = id;
	}
};
class CD8p : public CELL {
public:
	CD8p( PointF P = PointF( ) ) {
		VEC = PointF( );
		POS = P;
		DIV_INTERVAL = 0;
		MITOS = 0;
		LIFE = 20160;
		DIV_INTERVAL_CONST = 420;
		MITOS_CONST = 60;
		LIFE_CONST = 20160;
		SPEED = 10;
		RADIUS = 3.;
		ID = ID_CD8p;
	}
	bool can_move( ) {
		return MITOS == 0;
	}
};
class FLUID_ZONE {
public:
	FLUID_ZONE_TYPE type;
	CELL* myCell;
	double value [ FLUID_ZONE_N ];
	FLUID_ZONE( ) {
		type = NORMAL;
		for ( long i = 0; i < FLUID_ZONE_N; i++ ) {
			value [ i ] = 0;
		}
	}
};
struct FLUID_ZONE_DUP {
	double value [ FLUID_ZONE_N ];
};
PointL MAPSIZE = PointL( PointF( ) );
PointL MAP_MOVE = PointL( PointF( ) );
FLUID_ZONE**** MAP = NULL;
FLUID_ZONE_DUP**** dMAP = NULL;
map<CELL_ID , vector<bool> > semaphores_attract;
map<CELL_ID , vector<bool> > semaphores_fear;
map<CELL_ID , vector<double> > semaphores_produce;
void initSemaphores( ) {
	vector<bool> F;
	F.resize( FLUID_ZONE_N );
	// cd4+
	F [ GROWTH_FACTOR ] = true;
	semaphores_attract [ ID_CD4p ] = F;
	semaphores_attract [ ID_CD4pi ] = F;
	semaphores_attract [ ID_CD8p ] = F;
	F [ GROWTH_FACTOR ] = false;
	semaphores_fear [ ID_CD4p ] = F;
	semaphores_fear [ ID_CD4pi ] = F;
	semaphores_fear [ ID_CD8p ] = F;
	vector<double> F2;
	F2.resize( FLUID_ZONE_N );
	F2 [ GROWTH_FACTOR ] = 0.;
	semaphores_produce [ ID_CD4p ] = F2;
	semaphores_produce [ ID_CD4pi ] = F2;
	semaphores_produce [ ID_CD8p ] = F2;
}
PointL initMAP( vector<PointL> &P ) {
	long x = 0 , X = 0;
	long y = 0 , Y = 0;
	long z = 0 , Z = 0;
	if ( P.size( ) > 0 ) {
		x = X = P [ 0 ].x;
		y = Y = P [ 0 ].y;
		z = Z = P [ 0 ].z;
		for ( long i = 1; i < P.size( ); i++ ) {
			if ( P [ i ].x < x ) {
				x = P [ i ].x;
			}
			if ( P [ i ].y < y ) {
				y = P [ i ].y;
			}
			if ( P [ i ].z < z ) {
				z = P [ i ].z;
			}
			if ( P [ i ].x > X ) {
				X = P [ i ].x;
			}
			if ( P [ i ].y > Y ) {
				Y = P [ i ].y;
			}
			if ( P [ i ].z > Z ) {
				Z = P [ i ].z;
			}
		}
	}
	PointL R( -x , -y , -z );
	X -= x;
	Y -= y;
	Z -= z;
	X++;
	Y++;
	Z++;
	MAPSIZE = PointL( X , Y , Z );
	MAP = new FLUID_ZONE*** [ X ];
	dMAP = new FLUID_ZONE_DUP*** [ X ];
#pragma omp parallel for
	for ( long i = 0; i < X; i++ ) {
		MAP [ i ] = new FLUID_ZONE** [ Y ];
		dMAP [ i ] = new FLUID_ZONE_DUP** [ Y ];
		for ( long j = 0; j < Y; j++ ) {
			MAP [ i ] [ j ] = new FLUID_ZONE* [ Z ];
			dMAP [ i ] [ j ] = new FLUID_ZONE_DUP* [ Z ];
			for ( long k = 0; k < Z; k++ ) {
				MAP [ i ] [ j ] [ k ] = NULL;
				dMAP [ i ] [ j ] [ k ] = NULL;
			}
		}
	}
#pragma omp parallel for
	for ( long i = 0; i < P.size( ); i++ ) {
		MAP [ P [ i ].x - x ] [ P [ i ].y - y ] [ P [ i ].z - z ] = new FLUID_ZONE;
		dMAP [ P [ i ].x - x ] [ P [ i ].y - y ] [ P [ i ].z - z ] = new FLUID_ZONE_DUP;
	}
	return R;
}
bool check( PointL P ) {
	bool R = ( ( P.x > -1 ) && ( P.y > -1 ) && ( P.z > -1 ) );
	R &= P.x < MAPSIZE.x;
	R &= P.y < MAPSIZE.y;
	R &= P.z < MAPSIZE.z;
	if ( R ) {
		R = MAP [ P.x ] [ P.y ] [ P.z ] != NULL;
	}
	return R;
}
bool check( long x , long y , long z ) {
	bool R = ( ( x > -1 ) && ( y > -1 ) && ( z > -1 ) );
	R &= x < MAPSIZE.x;
	R &= y < MAPSIZE.y;
	R &= z < MAPSIZE.z;
	if ( R ) {
		R = MAP [ x ] [ y ] [ z ] != NULL;
	}
	return R;
}
bool checkfree( PointL P ) {
	bool R = check( P );
	if ( R ) {
		R = MAP [ P.x ] [ P.y ] [ P.z ]->myCell == NULL;
	}
	return R;
}
inline FLUID_ZONE *getFLink( long X , long Y , long Z ) {
	FLUID_ZONE *R = NULL;
	if ( check( X , Y , Z ) ) {
		R = MAP [ X ] [ Y ] [ Z ];
	}
	return R;
}
inline FLUID_ZONE *getFLink( PointL P ) {
	return getFLink( P.x , P.y , P.z );
}
void MAP_diffusion( ) {
	const double MAPS = 1.; // MAP block size
	double dc = 1.; // diffusion coefficient
#pragma omp parallel for
	for ( long i = 0; i < MAPSIZE.x; i++ ) {
		for ( long j = 0; j < MAPSIZE.y; j++ ) {
			for ( long k = 0; k < MAPSIZE.z; k++ ) {
				if ( check( i , j , k ) ) {
					for ( long f = 0; f < FLUID_ZONE_N; f++ ) {
						MAP [ i ] [ j ] [ k ]->value [ f ] -=
							dt *
							MAP [ i ] [ j ] [ k ]->value [ f ] *
							degradation [ f ];
					}
				}
			}
		}
	}
#pragma omp parallel for
	for ( long i = 0; i < MAPSIZE.x; i++ ) {
		for ( long j = 0; j < MAPSIZE.y; j++ ) {
			for ( long k = 0; k < MAPSIZE.z; k++ ) {
				if ( check( i , j , k ) ) {
					for ( long f = 0; f < FLUID_ZONE_N; f++ ) {
						dMAP [ i ] [ j ] [ k ]->value [ f ] = 0;
						double DX = 0 , DY = 0 , DZ = 0;
						{
							// d2v/dx2
							double from = MAP [ i ] [ j ] [ k ]->value [ f ];
							if ( check( i - 1 , j , k ) ) {
								from = MAP [ i - 1 ] [ j ] [ k ]->value [ f ];
							}
							double to = MAP [ i ] [ j ] [ k ]->value [ f ];
							if ( check( i + 1 , j , k ) ) {
								to = MAP [ i + 1 ] [ j ] [ k ]->value [ f ];
							}
							DX = ( to + from - 2. * MAP [ i ] [ j ] [ k ]->value [ f ] ) / ( MAPS * MAPS );
						}
						{
							// d2v/dy2
							double from = MAP [ i ] [ j ] [ k ]->value [ f ];
							if ( check( i , j - 1 , k ) ) {
								from = MAP [ i ] [ j - 1 ] [ k ]->value [ f ];
							}
							double to = MAP [ i ] [ j ] [ k ]->value [ f ];
							if ( check( i , j + 1 , k ) ) {
								to = MAP [ i ] [ j + 1 ] [ k ]->value [ f ];
							}
							DY = ( to + from - 2. * MAP [ i ] [ j ] [ k ]->value [ f ] ) / ( MAPS * MAPS );
						}
						{
							// d2v/dz2
							double from = MAP [ i ] [ j ] [ k ]->value [ f ];
							if ( check( i , j , k - 1 ) ) {
								from = MAP [ i ] [ j ] [ k - 1 ]->value [ f ];
							}
							double to = MAP [ i ] [ j ] [ k ]->value [ f ];
							if ( check( i , j , k + 1 ) ) {
								to = MAP [ i ] [ j ] [ k + 1 ]->value [ f ];
							}
							DZ = ( to + from - 2. * MAP [ i ] [ j ] [ k ]->value [ f ] ) / ( MAPS * MAPS );
						}
						dMAP [ i ] [ j ] [ k ]->value [ f ] = dc *( DX + DY + DZ ) * dt;
					}
				}
			}
		}
	}
#pragma omp parallel for
	for ( long i = 0; i < MAPSIZE.x; i++ ) {
		for ( long j = 0; j < MAPSIZE.y; j++ ) {
			for ( long k = 0; k < MAPSIZE.z; k++ ) {
				if ( check( i , j , k ) ) {
					for ( long f = 0; f < FLUID_ZONE_N; f++ ) {
						MAP [ i ] [ j ] [ k ]->value [ f ] += dMAP [ i ] [ j ] [ k ]->value [ f ];
					}
				}
			}
		}
	}
}
long applyBounds( ) {
	long R = 0;
	vector<PointL> BD;
	bool D1 = ( ( MAPSIZE.x == 1 ) && ( MAPSIZE.y == 1 ) ) || ( ( MAPSIZE.z == 1 ) && ( MAPSIZE.y == 1 ) ) || ( ( MAPSIZE.x == 1 ) && ( MAPSIZE.z == 1 ) );
	bool D2 = ( MAPSIZE.x == 1 ) || ( MAPSIZE.y == 1 ) || ( MAPSIZE.z == 1 );
	long CSUM = int( D1 ) * 2 + int( D2 && ( !D1 ) ) * 8 + int( ( !D1 ) && ( !D2 ) ) * 26;
#pragma omp parallel for
	for ( long i = 0; i < MAPSIZE.x; i++ ) {
		for ( long j = 0; j < MAPSIZE.y; j++ ) {
			for ( long k = 0; k < MAPSIZE.z; k++ ) {
				if ( check( i , j , k ) ) {
					int C = 0;
					for ( short a = -1; a < 2; a++ ) {
						for ( short b = -1; b < 2; b++ ) {
							for ( short c = -1; c < 2; c++ ) {
								if ( PointF( a , b , c ) != PointF( ) ) {
									C += int( check( i + a , j + b , k + c ) );
								}
							}
						}
					}
					if ( C != CSUM ) {
						MAP [ i ] [ j ] [ k ]->type = BOUND;
						R++;
					}
				}
			}
		}
	}
	return R;
}
void initZone( vector<PointL>& VL ) {
	// ЭТУ ФУНКЦИЮ МОЖНО И НУЖНО РЕДАКТИРОВАТЬ
	string R;
	if ( 1 ) {
		ifstream fin( "sphere.pobj" , ifstream::in );
		while ( !fin.eof( ) ) {
			fin >> R;
			if ( R == "vx" ) {
				long x , y , z;
				fin >> x >> y >> z;
				VL.push_back( PointL( x , y , z ) );
			}
		}
		fin.close( );
	}
	else {
		if ( 0 ) {
			for ( int i = 0; i < 51; i++ ) {
				for ( int j = 0; j < 51; j++ ) {
					VL.push_back( PointL( i , j , 0 ) );
				}
			}
		}
		else {
			if ( 0 ) {
				for ( int i = 0; i < 61; i++ ) {
					for ( int j = 0; j < 61; j++ ) {
						for ( int k = 0; k < 61; k++ ) {
							VL.push_back( PointL( i , j , k ) );
						}
					}
				}
			}
		}
	}
}
class CYCLE {
	PointL IP;
	long X , Y , Z , CT;
public:
	CYCLE( ) {
		IP.x = IP.y = IP.z = 0;
		X = Y = Z = 0;
		CT = 0;
	}
	void init( PointL P ) {
		IP = P;
		CT = IP.x * IP.y * IP.z;
		X = Y = Z = 0;
	}
	void reInit( ) {
		CT = IP.x * IP.y * IP.z;
		X = Y = Z = 0;
	}
	bool run( ) {
		return  ( CT != 0 );
	}
	PointL get( ) {
		PointL R;
		R.x = X;
		R.y = Y;
		R.z = Z;
		return R;
	}
	void next( ) {
		if ( CT ) {
			CT--;
			Z++;
			if ( Z >= IP.z ) {
				Y++;
				Z = 0;
				if ( Y >= IP.y ) {
					X++;
					Y = 0;
					if ( X >= IP.x ) {
						X = 0;
					}
				}
			}
		}
	}
	bool endOfX( ) {
		return ( X == ( IP.x - 1 ) ) && ( IP.x > 1 );
	}
	bool endOfY( ) {
		return ( Y == ( IP.y - 1 ) ) && ( IP.y > 1 );
	}
	bool endOfZ( ) {
		return ( Z == ( IP.z - 1 ) ) && ( IP.z > 1 );
	}
	bool end( ) {
		return endOfX( ) || endOfY( ) || endOfZ( );
	}
	void get( long &x , long& y , long& z ) {
		x = X;
		y = Y;
		z = Z;
	}
};
void applyModel( string S , PointL MOVE , FLUID_ZONE_TYPE FZType ) {
	vector<PointL> apply;
	string R;
	ifstream fin( S , ifstream::in );
	while ( !fin.eof( ) ) {
		fin >> R;
		if ( R == "vx" ) {
			long x , y , z;
			fin >> x >> y >> z;
			apply.push_back( PointL( x + MOVE.x , y + MOVE.y , z + MOVE.z ) );
		}
	}
	fin.close( );
#pragma omp parallel for
	for ( long i = 0; i < apply.size( ); i++ ) {
		PointL Q = apply [ i ];
		if ( check( Q ) ) {
			if ( MAP [ Q.x ] [ Q.y ] [ Q.z ]->type == NORMAL ) {
				MAP [ Q.x ] [ Q.y ] [ Q.z ]->type = FZType;
			}
		}
	}
}
// CD4+, cd4+i, CD8+
long initNums [ ] = { 1 , 0 , 0 };
vector<CELL*> all_cells;
vector<PointF> splitVector( PointF from , PointF to ) {
	double L = ( to - from ).length( );
	long ones = 2 + long( L );
	vector<PointF> R;
	R.push_back( from );
	if ( ones > 2 ) {
		PointF DP = ( to - from ) / ( ones - 1 );
		for ( long i = 1; i < ( ones - 1 ); i++ ) {
			R.push_back( from + ( DP * i ) );
		}
	}
	R.push_back( to );
	return R;
}
double try_move( CELL* cell ) {
	/*if ( cell->VEC.length( ) > EPS ) {
		cell->VEC = cell->VEC.normalize( );
	}*/
	cell->VEC -= cell->VEC * ( 1. - VEC_MOD ) * dt;
	PointF RMOD = PointF( rnd( ) , rnd( ) , rnd( ) ).normalize( ) * R_MOD * dt; // Random MODificator
	cell->VEC += RMOD;
	PointF move_to = cell->POS + cell->VEC * cell->SPEED * dt;
	/*cout << "Try move from (" << cell->POS.x << ", " << cell->POS.y << ", " << cell->POS.z << ") to (" <<
		move_to.x << ", " << move_to.y << ", " << move_to.z << ")" << endl;*/
	auto check_points = splitVector( cell->POS , move_to );
	double ML = 1.;
	for ( long i = 1; ( i < check_points.size( ) ) && ( ML == 1. ); i++ ) {
		PointL PL = check_points [ i ];
		if ( check( PL ) ) {
			bool Interruptor_Reasons = false;
			if ( getFLink( PL )->type != NORMAL ) {
				// врежемся в какую-нибудь структуру
				/*cout << "INTERRUPTED BY WALL" << endl;*/
				Interruptor_Reasons = true;
			}
			else {
				if ( getFLink( PL )->myCell != NULL && getFLink( PL )->myCell != cell ) {
					/*cout << "INTERRUPTED BY CELL" << endl;*/
					// врежемся в какую-то клетку
					Interruptor_Reasons = true;
				}
			}
			if ( Interruptor_Reasons ) {
				ML = double( i - 1 ) / double( check_points.size( ) - 1 );
			}
		}
		else {
			// вылетели из рабочей области
			ML = double( i - 1 ) / double( check_points.size( ) - 1 );
		}
	}
	if ( ML == 0 ) {
		PointF move_to_ [ 3 ] = { PointF( move_to.x , cell->POS.y , cell->POS.z ),
		PointF( cell->POS.x , move_to.y , cell->POS.z ),
		PointF( cell->POS.x , cell->POS.y , move_to.z ) };
		double ML_ [ 3 ];
		for ( int MLC = 0; MLC < 3; MLC++ ) {
			PointF move_to = move_to_ [ MLC ];
			/*cout << "Try move from (" << cell->POS.x << ", " << cell->POS.y << ", " << cell->POS.z << ") to (" <<
				move_to.x << ", " << move_to.y << ", " << move_to.z << ")" << endl;*/
			auto check_points = splitVector( cell->POS , move_to );
			double ML = 1.;
			for ( long i = 1; ( i < check_points.size( ) ) && ( ML == 1. ); i++ ) {
				PointL PL = check_points [ i ];
				if ( check( PL ) ) {
					bool Interruptor_Reasons = false;
					if ( getFLink( PL )->type != NORMAL ) {
						// врежемся в какую-нибудь структуру
						/*cout << "INTERRUPTED BY WALL" << endl;*/
						Interruptor_Reasons = true;
					}
					else {
						if ( getFLink( PL )->myCell != NULL && getFLink( PL )->myCell != cell ) {
							/*cout << "INTERRUPTED BY CELL" << endl;*/
							// врежемся в какую-то клетку
							Interruptor_Reasons = true;
						}
					}
					if ( Interruptor_Reasons ) {
						ML = double( i - 1 ) / double( check_points.size( ) - 1 );
					}
				}
				else {
					// вылетели из рабочей области
					ML = double( i - 1 ) / double( check_points.size( ) - 1 );
				}
			}
			ML_ [ MLC ] = ML;
		}
		double MLL = ML;
		PointF move_too = move_to;
		for ( int MLC = 0; MLC < 3; MLC++ ) {
			if ( ( ( move_to_ [ MLC ] - cell->POS ) * ML_ [ MLC ] ) > ( ( move_too - cell->POS ) * MLL ) ) {
				MLL = ML_ [ MLC ];
				move_too = move_to_ [ MLC ];
			}
		}
		ML = MLL;
		move_to = move_too;
	}
	PointF final_pos = cell->POS + ( move_to - cell->POS ) * ML;
	// мы провели все необходимые проверки, теперь можно сменить координату клетки
	double R = ( cell->POS - final_pos ).length( );
	PointL WAS_FLUID_ZONE( cell->POS );
	PointL BECAME_FLUID_ZONE( final_pos );
	if ( !cell->MITOS ) {
		getFLink( WAS_FLUID_ZONE )->myCell = NULL;
		getFLink( BECAME_FLUID_ZONE )->myCell = cell;
		cell->POS = final_pos;
	}
	return R;
}
void throw_fluid( PointL p , long ID , double volume ) {
	if ( check( p ) && ID >= 0 && ID < FLUID_ZONE_N && volume >= 0. ) {
		MAP [ p.x ] [ p.y ] [ p.z ]->value [ ID ] += volume;
	}
}
void throw_fluid( CELL* cell , long ID , double volume ) {
	if ( cell != NULL && ID >= 0 && ID < FLUID_ZONE_N && volume >= 0. ) {
		vector<PointL> to_throw;
		PointL CC( cell->POS );
		long CR = long( cell->RADIUS ) + 1;
		for ( long dx = -CR; dx <= CR; dx++ ) {
			long L1 = long( pow( CR * CR - dx * dx , 0.5 ) ) + 1;
			for ( long dy = -L1; dy <= L1; dy++ ) {
				long L2 = int( ( CR * CR - dx * dx - dy * dy ) > 0 );
				if ( L2 ) {
					L2 = long( pow( CR * CR - dx * dx - dy * dy , 0.5 ) ) + 1;
					for ( long dz = -L2; dz <= L2; dz++ ) {
						PointL CPL( CC.x + dx , CC.y + dy , CC.z + dz );
						double dist = ( cell->POS - PointF( CC.x + dx , CC.y + dy , CC.z + dz ) ).length( );
						if ( dist <= cell->RADIUS && check( CPL ) ) {
							to_throw.push_back( CPL );
						}
					}
				}
			}
		}
		if ( to_throw.size( ) ) {
			// размажем произведённое вещество по всему пространству, занимаемому клеткой
			volume /= to_throw.size( );
			for ( long i = 0; i < to_throw.size( ); i++ ) {
				MAP [ to_throw [ i ].x ] [ to_throw [ i ].y ] [ to_throw [ i ].z ]->value [ ID ] += volume;
			}
		}
	}
}
vector<PointL> grabZonesOfCell( CELL* cell ) {
	vector<PointL> to_analyse;
	PointL CC( cell->POS );
	long CR = long( MAX_CELL_DIAMETER ) + 1;
	for ( long dx = -CR; dx <= CR; dx++ ) {
		long L1 = long( pow( CR * CR - dx * dx , 0.5 ) ) + 1;
		for ( long dy = -L1; dy <= L1; dy++ ) {
			long L2 = int( ( CR * CR - dx * dx - dy * dy ) > 0 );
			if ( L2 ) {
				L2 = long( pow( CR * CR - dx * dx - dy * dy , 0.5 ) ) + 1;
				for ( long dz = -L2; dz <= L2; dz++ ) {
					//getchar( );
					PointL CPL( CC.x + dx , CC.y + dy , CC.z + dz );
					double dist = ( cell->POS - CPL.convert( ) ).length( );
					if ( check( CPL ) ) {
						if ( getFLink( CPL )->type == NORMAL ) {
							to_analyse.push_back( CPL );
						}
					}
				}
			}
		}
	}
	return to_analyse;
}
PointF selectMovePath( CELL* cell , long j/*fluid id*/ ) {
	vector<PointL> to_analyse = grabZonesOfCell( cell );
	vector<double> FLUIDS;
	long FL_MIN = 0;
	double valMax = 0;
	long VMAX = 0;
	for ( long i = 0; i < to_analyse.size( ); i++ ) {
		PointL PA = to_analyse [ i ];
		FLUIDS.push_back( getFLink( PA )->value [ j ] );
		if ( valMax < FLUIDS [ i ] ) {
			valMax = FLUIDS [ i ];
			VMAX = i;
		}
		if ( getFLink( PA )->value [ j ] < getFLink( to_analyse [ FL_MIN ] )->value [ j ] ) {
			FL_MIN = i;
		}
	}
	cout << "VALMAX = " << valMax << endl;
	double SUM = 0.;
	double minus_all = FLUIDS [ FL_MIN ];
	cout << "MINUS_ALL = " << minus_all << endl;
	for ( long i = 0; i < FLUIDS.size( ); i++ ) {
		FLUIDS [ i ] -= minus_all;
		SUM += FLUIDS [ i ];
	}
	cout << "SUM = " << SUM << endl;
	SUM *= RND( );
	PointF NEW_PATH( cell->POS );
	if ( SUM > EPS ) {
		long PATHGO = -1;
		while ( SUM > 0. ) {
			PATHGO++;
			SUM -= FLUIDS [ PATHGO ];
		}
		NEW_PATH = to_analyse [ 0 * VMAX + 1 * PATHGO ].convert( );
	}
	return NEW_PATH;
}
double grabFluidIntencity( CELL* cell , long j/*fluid id*/ ) {
	vector<PointL> to_analyse = grabZonesOfCell( cell );
	vector<double> grabbed;
	double sum = 0.;
	for ( long i = 0; i < to_analyse.size( ); i++ ) {
		grabbed.push_back( getFLink( to_analyse [ i ] )->value [ j ] );
		sum += grabbed [ i ];
	}
	if ( grabbed.size( ) > 1 ) {
		sum /= grabbed.size( );
	}
	return sum;
}
PointL pex1 , pex2;
void cells_dynamic( ) {
#pragma omp parallel for
	for ( long i = 0; i < all_cells.size( ); i++ ) {
		srand( unsigned int( i ) );
		bool end_of_mitos;
		// обновляем возраст
		all_cells [ i ]->stat_time_update( dt , end_of_mitos );
		if ( all_cells [ i ]->live( ) ) {
			if ( all_cells [ i ]->ID == ID_CD4p || all_cells [ i ]->ID == ID_CD4pi ) {
				// делаем указатель на типизированную клетку
				CD4p *cell = ( CD4p* ) all_cells [ i ];
				if ( end_of_mitos ) {
					// закончился митоз, создаём клетку
					CELL* new_cell = new CD4p( cell->POS + noize( ) , cell->ID );
#pragma omp critical
					all_cells.push_back( new_cell );
				}
			}
			if ( all_cells [ i ]->ID == ID_CD8p ) {
				// делаем указатель на типизированную клетку
				CD8p *cell = ( CD8p* ) all_cells [ i ];
				if ( end_of_mitos ) {
					// закончился митоз, создаём клетку
					CELL* new_cell = new CD8p( cell->POS + noize( ) );
#pragma omp critical
					all_cells.push_back( new_cell );
				}
			}
			if ( end_of_mitos ) {
				// восстановить время жизни у старой клетки и задержку деления
				all_cells [ i ]->DIV_INTERVAL = all_cells [ i ]->DIV_INTERVAL_CONST;
				all_cells [ i ]->LIFE = all_cells [ i ]->LIFE_CONST;
			}
			try_move( all_cells [ i ] );
			/*cout <<
				"CELL[" << i << "], " <<
				"POS=(" << all_cells [ i ]->POS.x << ", " << all_cells [ i ]->POS.y << ", " << all_cells [ i ]->POS.z << "), " <<
				"VEC=(" << all_cells [ i ]->VEC.x << ", " << all_cells [ i ]->VEC.y << ", " << all_cells [ i ]->VEC.z << ") [" <<
				all_cells [ i ]->VEC.length( ) << "]" << endl;*/
		}
	}
	for ( long i = 0; i < all_cells.size( ); i++ ) {
		if ( all_cells [ i ]->live( ) ) {
			// сканируем пространство вокруг
			CELL *cell = all_cells [ i ];
			vector<CELL*> to_interact;
			vector<PointL> to_push_away;
			PointL CC( cell->POS );
			long CR = long( MAX_CELL_DIAMETER ) + 1;
			// поиск препятствий, отталкивание от них, взаимодействие с клетками
			for ( long dx = -CR; dx <= CR; dx++ ) {
				long L1 = long( pow( CR * CR - dx * dx , 0.5 ) ) + 1;
				for ( long dy = -L1; dy <= L1; dy++ ) {
					long L2 = int( ( CR * CR - dx * dx - dy * dy ) > 0 );
					if ( L2 ) {
						L2 = long( pow( CR * CR - dx * dx - dy * dy , 0.5 ) ) + 1;
						for ( long dz = -L2; dz <= L2; dz++ ) {
							PointL CPL( CC.x + dx , CC.y + dy , CC.z + dz );
							double dist = ( cell->POS - CPL.convert( ) ).length( );
							if ( check( CPL ) ) {
								if ( getFLink( CPL )->type != NORMAL || getFLink( CPL )->type != VESSEL ) {
									to_push_away.push_back( CPL );
								}
								if ( getFLink( CPL )->myCell != NULL ) {
									if ( getFLink( CPL )->myCell->live( ) && getFLink( CPL )->myCell != cell ) {
										to_interact.push_back( getFLink( CPL )->myCell );
									}
								}
							}
						}
					}
				}
			}
			PointF PushAwayVEC = PointF( );
			for ( long j = 0; j < to_push_away.size( ); j++ ) {
				PointF TOUCH( to_push_away [ j ].x , to_push_away [ j ].y , to_push_away [ j ].z );
				double D = ( cell->POS - TOUCH ).length( );
				if ( D < cell->RADIUS ) {
					// вообще-то, мы скребём стенку боком
					if ( D < EPS ) {
						D = EPS;
						PointF goAway = PointF( cell->POS - TOUCH ).normalize( ) * ( cell->RADIUS - D );
						PushAwayVEC += goAway;
					}
				}
			}
			if ( PushAwayVEC.length( ) > 1. ) {
				PushAwayVEC = PushAwayVEC.normalize( );
			}
			cell->VEC += PushAwayVEC * WALL__TOUCH_FORCE_COEFFICIENT * dt;
			PointF InteractVEC = PointF( );
			for ( long j = 0; j < to_interact.size( ); j++ ) {
				CELL* cell2 = to_interact [ j ];
				double D = ( cell->POS - cell2->POS ).length( );
				if ( D <= ( cell->RADIUS + cell2->RADIUS ) ) {
					if ( D < EPS ) {
						D = EPS;
					}
					PointF goAway = PointF( cell->POS - cell2->POS ).normalize( ) * ( cell->RADIUS + cell2->RADIUS - D ) * 0.5/*та клетка тоже оттолкнётся*/;
					InteractVEC += goAway;
					{
						/*
							ЗДЕСЬ ДОБАВЬТЕ ВСЯКУЮ ФИГНЮ ТИПА УБИЙСТВА КЛЕТОК И Т.П.
						*/
					}
				}
			}
			if ( InteractVEC.length( ) > 1. ) {
				InteractVEC = InteractVEC.normalize( );
			}
			cell->VEC += InteractVEC * CELL__TOUCH_FORCE_COEFFICIENT * dt;
			// взаимодействие с флюидами
			PointF FluidVEC = PointF( );
			for ( long j = 0; j < FLUID_ZONE_N; j++ ) {
				PointF attr_point = selectMovePath( cell , j );
				PointF attr_vec = attr_point - cell->POS;
				attr_vec *= int( semaphores_attract [ cell->ID ] [ j ] ) - int( semaphores_fear [ cell->ID ] [ j ] );
				FluidVEC += attr_vec.normalize( );
				if ( j == GROWTH_FACTOR ) {
					double grab = grabFluidIntencity( cell , j );
					if ( cell->ID == ID_CD4p || cell->ID == ID_CD4pi ) {
						if ( grab > 60 ) { // 60%
							if ( cell->can_div( ) ) {
								cell->MITOS = cell->MITOS_CONST;
							}
						}
					}
				}
			}
			if ( FluidVEC.length( ) > 1. ) {
				FluidVEC = FluidVEC.normalize( );
			}
			cell->VEC += FluidVEC * FLUID_TOUCH_FORCE_COEFFICIENT * dt;
			cout << "CELL POS = (" << cell->POS.x << ", " << cell->POS.y << ", " << cell->POS.z << ")" << endl;
			cout << "VEC TO SOURCE = (" << pex1.convert( ).normalize( ).x << ", " << pex1.convert( ).normalize( ).y << ", " << pex1.convert( ).normalize( ).z << ")" << endl;
			cout << "COS of ANGLES = " << ( pex1.convert( ).normalize( ) | FluidVEC.normalize( ) ) << endl;
			//getchar( );
		}
	}
	// удаление мёртвых клеток
#pragma omp parallel for
	for ( long i = 0; i < all_cells.size( ); i++ ) {
		if ( all_cells [ i ]->dead( ) ) {
			if ( all_cells [ i ]->ID == ID_CD4p ) {
				delete ( CD4p* ) all_cells [ i ];
				all_cells [ i ] = NULL;
			}
			else {
				if ( all_cells [ i ]->ID == ID_CD8p ) {
					delete ( CD8p* ) all_cells [ i ];
					all_cells [ i ] = NULL;
				}
			}
		}
	}
	// удаляем нуль-указатели
	for ( long i = long( all_cells.size( ) ) - 1; i >= 0; i-- ) {
		if ( all_cells [ i ] == NULL ) {
			all_cells [ i ] = all_cells [ all_cells.size( ) ];
			all_cells.pop_back( );
		}
	}
}

template <typename T> void placeCellInRandomPlace( T ) {
	CELL* cell = new T( );
	bool canPlace = false;
	long a , b , c;
	while ( !canPlace ) {
		a = rand( ) % MAPSIZE.x;
		b = rand( ) % MAPSIZE.y;
		c = rand( ) % MAPSIZE.z;
		if ( check( a , b , c ) ) {
			if ( checkfree( PointL( a , b , c ) ) ) {
				canPlace = true;
			}
		}
	}
	all_cells.push_back( cell );
	cell->POS = PointF( a , b , c );
	getFLink( a , b , c )->myCell = cell;
}
template <typename T> bool placeCellInPlace( T , PointL P ) {
	bool canPlace = false;
	long a , b , c;
	a = rand( ) % MAPSIZE.x;
	b = rand( ) % MAPSIZE.y;
	c = rand( ) % MAPSIZE.z;
	a = P.x;
	b = P.y;
	c = P.z;
	if ( check( a , b , c ) ) {
		CELL* cell = new T( );
		if ( checkfree( PointL( a , b , c ) ) ) {
			canPlace = true;
		}
		if ( canPlace ) {
			all_cells.push_back( cell );
			cell->POS = PointF( a , b , c );
			getFLink( a , b , c )->myCell = cell;
		}
		else {
			delete cell;
		}
	}
	return canPlace;
}

PointF temperatureColor( double value , double max , double min ) {
	PointF R = PointF( );
	if ( ( max - min ) < EPS ) {
		return R;
	}
	double D = max - min;
	PointF MAX( 1 , 0 , 0 );
	PointF MID( 0 , 1 , 0 );
	PointF MIN( 0 , 0 , 1 );
	if ( value > max ) {
		value = max;
	}
	if ( value < min ) {
		value = min;
	}
	if ( value >= ( max - D / 2. ) ) {
		R = MAX * ( 1. - ( max - value ) / ( D / 2. ) ) + MID * ( max - value ) / ( D / 2. );
	}
	else {
		max = max - D / 2.;
		R = MID * ( 1. - ( max - value ) / ( D / 2. ) ) + MIN * ( max - value ) / ( D / 2. );
	}
	return R;
}

double anaFunc( double x , double t ) {
	double sum = 0.;
	double psum = 0.;
	long one = -1;
	while ( one == -1 || psum > EPS || one < 10000 ) {
		one++;
		psum = 4 / PI * exp( -PI * PI * ( 2 * one + 1 )*( 2 * one + 1 ) * t ) / ( 2 * one + 1 );
		sum += psum * sin( PI * ( 2 * one + 1 ) * x );
	}
	return sum;
}
double ana3Dfunc( double x , double y , double z , double t , bool cou = false ) {
	if ( cou ) {
		cout << "X = " << x << endl;
		cout << "Y = " << y << endl;
		cout << "Z = " << z << endl;
		cout << "T = " << t << endl;
		cout << "First component (1 / (4*pi*t)^(3/2)) = " << ( 0.125 / ( pow( 4. * PI * t , 1.5 ) ) ) << endl;
		cout << "Second component (e^(-(1 / 4t) * (x^2+y^2+z^2))) = " << exp( ( -1. / ( 4. * t ) ) * ( x * x + y * y + z * z ) ) << endl;
	}
	return ( 0.125 / ( pow( 4. * PI * t , 1.5 ) ) )*exp( ( -1. / ( 4. * t ) ) * ( x * x + y * y + z * z ) );
}
int main( int argc , char** argv ) {
	lmModel lmcube , lmcell;
	lmcube.loadModel( "mod_cube.obj" );
	lmcell.loadModel( "mod_cell.obj" );
	if ( 1 ) {
		// вся рабочая область
		//system( "vpc.exe -ss 1.0 -lo sphere10R.obj -ff -a -cl -wq sphere.pobj" );
		// вся рабочая область
		//system( "vpc.exe -ss 1.0 -lo sphere50RXXS.obj -ff -a -cl -wq sphere.pobj" );
		// часть, которую нужно отпилить от ФРК
		//system( "vpc.exe -ss 1.0 -lo tffvpss.obj -ff -a -cl -lv sphere.pobj -r -cl -wq remove.pobj" );
		// часть, которую нужно использовать в сети ФРК
		//system( "vpc.exe -ss 1.0 -lo tffvpss.obj -ff -a -cl -lv remove.pobj -r -cl -wq frc_part.pobj" );
	}
	{
		// НАСТРОЙТЕ ЗДЕСЬ КОЭФФИЦИЕНТЫ ДЕГРАДАЦИИ ВЕЩЕСТВ
		degradation [ GROWTH_FACTOR ] = 0.2;
	}
	vector<PointL> VL;
	initZone( VL );
	initSemaphores( );
	MAP_MOVE = initMAP( VL );
	long count = 0;
	applyBounds( );
	applyModel( "frc_part.pobj" , MAP_MOVE , FRC );
	vector<PointL> toPushFRC , toPushVESSEL;
	for ( long i = 0; i < MAPSIZE.x; i++ ) {
		for ( long j = 0; j < MAPSIZE.y; j++ ) {
			for ( long k = 0; k < MAPSIZE.z; k++ ) {
				if ( check( i , j , k ) ) {
					if ( getFLink( i , j , k )->type == FRC ) {
						toPushFRC.push_back( PointL( i , j , k ) );
					}
					if ( getFLink( i , j , k )->type == VESSEL ) {
						toPushVESSEL.push_back( PointL( i , j , k ) );
					}
				}
			}
		}
	}
	{
		// внесение клеток в модель
		for ( int i = 0; i < 700; i++ ) {
			placeCellInRandomPlace( CD4p( ) );
		}
		/*
			Разработка агентной модели миграции и взаимодействия клеточных популяций и ВИЧ в замкнутой области лимфатического узла
		*/
	}
	gtc_start( );
	double timeInfo = omp_get_wtime( ) , fullTime = 0;
	system( "del cells_pos_view*.obj" );
	system( "del screen*.bmp" );
	for ( long SRETI = 0; SRETI < ITERS; SRETI++ ) {
		{
			// установление граничных условий
#pragma omp parallel for
			for ( long i = 0; i < toPushFRC.size( ); i++ ) {
				//getFLink( toPushFRC [ i ] )->value [ GROWTH_FACTOR ] = 30. + 50. * int( toPushFRC [ i ].x == 0 );
				getFLink( toPushFRC [ i ] )->value [ GROWTH_FACTOR ] = 100.;
			}
#pragma omp parallel for
			for ( long i = 0; i < toPushVESSEL.size( ); i++ ) {
				getFLink( toPushVESSEL [ i ] )->value [ GROWTH_FACTOR ] = 0.;
			}
			// диффузия в среде
			MAP_diffusion( );
			// установление граничных условий
#pragma omp parallel for
			for ( long i = 0; i < toPushFRC.size( ); i++ ) {
				//getFLink( toPushFRC [ i ] )->value [ GROWTH_FACTOR ] = 30. + 50. * int( toPushFRC [ i ].x == 0 );
				getFLink( toPushFRC [ i ] )->value [ GROWTH_FACTOR ] = 100.;
			}
#pragma omp parallel for
			for ( long i = 0; i < toPushVESSEL.size( ); i++ ) {
				getFLink( toPushVESSEL [ i ] )->value [ GROWTH_FACTOR ] = 0.;
			}
		}
		// динамика клеток
		cells_dynamic( );
		cout << endl << endl;
		// Захват времени выполнения, вычисление оставшегося времени
		if ( omp_get_wtime( ) - timeInfo > 10. ) {
			timeInfo = omp_get_wtime( ) - timeInfo;
			double perIter = ( fullTime + timeInfo ) / ( SRETI + 1 );
			cout << "Please wait ";
			double waitSec = perIter * ( ITERS - ( SRETI + 1 ) );
			long waitMin = long( waitSec / 60 );
			long waitHours = waitMin / 60;
			long waitDays = waitHours / 24;
			long waitYears = long( waitDays / 365 );
			if ( waitYears ) {
				cout << waitYears << " years, ";
				waitDays = waitDays % 365;
			}
			if ( waitDays ) {
				cout << waitDays << " days, ";
				waitHours = waitHours % 24;
			}
			if ( waitHours ) {
				cout << waitHours << " hours, ";
				waitMin = waitMin % 60;
			}
			if ( waitMin ) {
				cout << waitMin << " minutes, ";
				waitSec = waitSec - long( waitSec / 60 ) * 60;
			}
			cout << long( waitSec ) << " seconds";
			cout << " ( " << fullTime + timeInfo << " sec spent )" << endl;
			fullTime += timeInfo;
			timeInfo = omp_get_wtime( );
		}
		bool REPORTS = true;
		if ( REPORTS ) {
			double DX = -MAPSIZE.x / 2.;
			double DY = -MAPSIZE.y / 2.;
			double DZ = -MAPSIZE.z / 2.;
			PointF DXYZ( DX , DY , DZ );
			char WART [ 4096 ];
			sprintf( WART , "_%07ld" , SRETI );
			string WARTS( WART );
			lmpoints.clear( );
			lmtries.clear( );
			lmColors.clear( );
			long ct = 0;
			for ( long i = 0; i < MAPSIZE.x; i++ ) {
				for ( long j = 0; j < MAPSIZE.y; j++ ) {
					for ( long k = 0; k < MAPSIZE.z; k++ ) {
						if ( check( i , j , k ) ) {
							if ( getFLink( i , j , k )->type != NORMAL &&  getFLink( i , j , k )->type != BOUND ) {
								/*loadModelForRezult( "mod_cube.obj" ,
													PointF( i , j , k ) + DXYZ ,
													0.125 * ( 1 + int( getFLink( i , j , k )->type == FRC ) * 7 ) ,
													PointF(
													int( getFLink( i , j , k )->type != FRC ) ,
													1 ,
													int( getFLink( i , j , k )->type != FRC )
								)
								);*/
								pushLmModel( lmcube , PointF( i , j , k ) + DXYZ ,
											 0.125 * ( 1 + int( getFLink( i , j , k )->type == FRC ) * 7 ) ,
											 PointF(
											 int( getFLink( i , j , k )->type != FRC ) ,
											 1 ,
											 int( getFLink( i , j , k )->type != FRC )
								) * 0.5 );
							}
							if ( false && getFLink( i , j , k )->type == NORMAL &&  getFLink( i , j , k )->value [ GROWTH_FACTOR ] > EPS ) {
								/*loadModelForRezult(
									"mod_cube.obj" ,
									PointF( i , j , k ) + DXYZ ,
									0.0099 * getFLink( i , j , k )->value [ GROWTH_FACTOR ] ,
									temperatureColor( getFLink( i , j , k )->value [ GROWTH_FACTOR ] , 100. , 0. )
								);*/
								pushLmModel(
									lmcube ,
									PointF( i , j , k ) + DXYZ ,
									0.0099 * getFLink( i , j , k )->value [ GROWTH_FACTOR ] ,
									temperatureColor( getFLink( i , j , k )->value [ GROWTH_FACTOR ] , 100. , 0. )
								);
							}
						}
					}
				}
			}
			for ( long i = 0; i < all_cells.size( ); i++ ) {
				//loadModelForRezult( "mod_cell.obj" , all_cells [ i ]->POS + DXYZ , all_cells [ i ]->RADIUS , PointF( int( all_cells [ i ]->ID == ID_CD8p ) , int( all_cells [ i ]->ID == ID_CD4p ) , int( all_cells [ i ]->ID == ID_CD4pi ) ) );
				pushLmModel( lmcell , all_cells [ i ]->POS + DXYZ , all_cells [ i ]->RADIUS , PointF( int( all_cells [ i ]->ID == ID_CD8p ) , int( all_cells [ i ]->ID == ID_CD4p ) , int( all_cells [ i ]->ID == ID_CD4pi ) ) );
			}
			throwAllModels2Out( ( "cells_pos_view" + WARTS + ".obj" ).c_str( ) );
			system( ( "OPOvis.exe -i cells_pos_view" + WARTS + ".obj -iq FRCnvi.obj -noslice -o screen" + WARTS + ".bmp -md 100.0 -rx 0.0 -ry " + "90.0" + " -rz 0.0" ).c_str( ) );
		}
		//getchar( );
	}
	system( "del video.mp4" );
	system( "ffmpeg.exe -r 10/1 -i screen_%07d.bmp -c:v libx264 -vf \"fps = 50 , format = yuv420p\" video.mp4" );
	system( "del cells_pos_view*.obj" );
	system( "del screen*.bmp" );
	timeInfo = gtc_stop( );
	cout << "Time spent " << timeInfo << " seconds." << endl;
	cout << timeInfo / ( ITERS + int( ITERS == 0 ) ) << " sec per iteration" << endl;
	{
		ofstream foutX( "fluidX.txt" , ofstream::out );
		for ( long j = 0; j < MAPSIZE.y; j++ ) {
			for ( long i = 0; i < MAPSIZE.z; i++ ) {
				long k = 10;
				//long k = 0;//long( MAPSIZE.z / 2 );
				if ( check( k , j , i ) ) {
					foutX << /*anaFunc( j * 0.1 , 0.1 ) - */getFLink( k , j , i )->value [ GROWTH_FACTOR ];
				}
				else {
					foutX << "0.0";
				}
				if ( i != ( MAPSIZE.z - 1 ) ) {
					foutX << " ";
				}
			}
			foutX << endl;
		}
		foutX.close( );
		ofstream foutY( "fluidY.txt" , ofstream::out );
		for ( long j = 0; j < MAPSIZE.x; j++ ) {
			for ( long i = 0; i < MAPSIZE.z; i++ ) {
				long k = 10;
				//long k = 0;//long( MAPSIZE.z / 2 );
				if ( check( j , k , i ) ) {
					foutY << /*anaFunc( j * 0.1 , 0.1 ) - */getFLink( j , k , i )->value [ GROWTH_FACTOR ];
				}
				else {
					foutY << "0.0";
				}
				if ( i != ( MAPSIZE.z - 1 ) ) {
					foutY << " ";
				}
			}
			foutY << endl;
		}
		foutY.close( );
		ofstream foutZ( "fluidZ.txt" , ofstream::out );
		for ( long j = 0; j < MAPSIZE.x; j++ ) {
			for ( long i = 0; i < MAPSIZE.y; i++ ) {
				long k = 30;
				//long k = 0;//long( MAPSIZE.z / 2 );
				if ( check( j , i , k ) ) {
					//foutZ << getFLink( j , i , k )->value [ GROWTH_FACTOR ];
					foutZ << ana3Dfunc( j * 0.5 - 15. , i * 0.5 - 15. , 0. , ft );
					//foutZ << ana3Dfunc( j * 0.5 - 15. , i * 0.5 - 15. , 0. , ft ) - getFLink( j , i , k )->value [ GROWTH_FACTOR ];
				}
				else {
					foutZ << "0.0";
				}
				if ( i != ( MAPSIZE.y - 1 ) ) {
					foutZ << " ";
				}
			}
			foutZ << endl;
		}
		foutZ.close( );
		getchar( );
	}
	// ### OUT DATA WRITE ###
	ofstream fout( "gnuplot_script" , ofstream::out );
	fout << "set term png size 1920,1080 enhanced font 'Verdana,24'" << endl << "set output \"fluidX.png\"" << endl << "set grid" << endl << "plot \"fluidX.txt\" matrix with image" << endl;// using 1 with lines;
	fout << "set term png size 1920,1080 enhanced font 'Verdana,24'" << endl << "set output \"fluidY.png\"" << endl << "set grid" << endl << "plot \"fluidY.txt\" matrix with image" << endl;// using 1 with lines;
	fout << "set term png size 1920,1080 enhanced font 'Verdana,24'" << endl << "set output \"fluidZ.png\"" << endl << "set grid" << endl << "plot \"fluidZ.txt\" matrix with image" << endl;// using 1 with lines;
	fout.close( );
	system( "gnuplot.exe gnuplot_script" );
	cout << MAPSIZE.x << " " << MAPSIZE.y << " " << MAPSIZE.z << endl;
	while ( getchar( ) != '\n' );
	return 0;
}