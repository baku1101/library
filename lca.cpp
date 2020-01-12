/*
The implementation of <O(n), O(1)> LCA with C++
varified with GRL_5_C

http://joisino.hatenablog.com/entry/2017/08/13/210000

Copyright (c) 2017 joisino
Released under the MIT license
http://opensource.org/licenses/mit-license.php
*/

template<class T>
class MinOp{
	public:
		T operator () ( T a , T b ){ return min( a , b ); }
};

// sparse table
// static range semigroup query
// time complexity: <O(n log n), O(1)>
// OpFunc is binary operator: T x T -> T
template<typename T, typename OpFunc>
struct SparseTable{
	OpFunc op;
	int size;
	vector<int> lg;
	vector<vector<T> > table;
	void init( const vector<T> &array, OpFunc opfunc ){
		int n = array.size();
		op = opfunc;

		lg.assign( n + 1 , 0 );
		for( int i = 1; i <= n; i++ ){
			lg[i] = 31 - __builtin_clz( i );
		}

		table.assign( lg[n] + 1, array );
		for( int i = 1; i <= lg[n]; i++ ){
			for( int j = 0; j < n; j++ ){
				if( j + (1<<(i-1)) < n ){
					table[i][j] = op( table[i-1][j] , table[i-1][j+(1<<(i-1))] );
				} else {
					table[i][j] = table[i-1][j];
				}
			}
		}
	}
	T query( int l , int r ){
		assert( l < r );
		return op( table[lg[r-l]][l], table[lg[r-l]][r-(1<<lg[r-l])] );
	}
};


// plus minus one Range Minimum Query
// time complexity: <O(n), O(1)>
struct PMORMQ{
	vector<int> a;
	SparseTable<pair<int,int>,MinOp<pair<int,int> > > sparse_table;
	vector<vector<vector<int> > > lookup_table;
	vector<int> block_type;
	int block_size, n_block;
	void init( const vector<int> &array ){
		a = array;
		int n = a.size();
		block_size = max( 1, ( 31 - __builtin_clz( n ) ) / 2 );
		while( n % block_size != 0 ){
			a.push_back( a.back() + 1 );
			n++;
		}
		n_block = n / block_size;

		vector<pair<int,int> > b( n_block, make_pair( INT_MAX, INT_MAX ) );
		for( int i = 0; i < n; i++ ){
			b[i/block_size] = min( b[i/block_size], make_pair( a[i], i ) );
		}
		sparse_table.init( b, MinOp<pair<int,int> >() );

		block_type.assign( n_block, 0 );
		for( int i = 0; i < n_block; i++ ){
			int cur = 0;
			for( int j = 0; j < block_size-1; j++ ){
				int ind = i * block_size + j;
				if( a[ind] < a[ind+1] ){
					cur |= 1 << j;
				}
			}
			block_type[i] = cur;
		}

		lookup_table.assign( 1 << (block_size-1), vector<vector<int> >( block_size, vector<int>( block_size+1 ) ) );
		for( int i = 0; i < (1<<(block_size-1)); i++ ){
			for( int j = 0; j < block_size; j++ ){
				int res = 0;
				int cur = 0;
				int pos = j;
				for( int k = j+1; k <= block_size; k++ ){
					lookup_table[i][j][k] = pos;
					if( i & ( 1 << (k-1) ) ){
						cur++;
					} else {
						cur--;
					}
					if( res > cur ){
						pos = k;
						res = cur;
					}
				}
			}
		}
	}
	int query( int l, int r ){ // return position
		assert( l < r );
		int lb = l / block_size;
		int rb = r / block_size;
		if( lb == rb ){
			return lb * block_size + lookup_table[ block_type[lb] ][ l % block_size ][ r % block_size ];
		}
		int pl = lb * block_size + lookup_table[ block_type[lb] ][ l % block_size ][ block_size ];
		int pr = rb * block_size + lookup_table[ block_type[rb] ][0][ r % block_size ];
		int pos = pl;
		if( r % block_size > 0 && a[pl] > a[pr] ){
			pos = pr;
		}
		if( lb + 1 == rb ){
			return pos;
		}
		int sp = sparse_table.query( lb+1, rb ).second;
		if( a[pos] > a[sp] ){
			return sp;
		}
		return pos;
	}
};

// LCA
// time complexity: <O(n), O(1)>
struct LCA{
	int n;
	vector<vector<int> > G;
	PMORMQ rmq;
	int cnt;
	vector<int> depth, id, in;
	void init( int size ){
		n = size;
		G.assign( n, vector<int>( 0 ) );
	}
	void add_edge( int a , int b ){
		G[a].push_back( b );
		G[b].push_back( a );
	}
	void dfs( int x , int p , int d ){
		id[cnt] = x;
		depth.push_back( d );
		in[x] = cnt++;
		for( int v : G[x] ){
			if( v == p ){
				continue;
			}
			dfs( v , x , d+1 );
			id[cnt] = x;
			depth.push_back( d );
			cnt++;
		}
	}
	void precalc( int root ){
		cnt = 0;
		depth.clear();
		in.assign( n, -1 );
		id.assign( 2*n-1, -1 );
		dfs( root, -1, 0 );
		rmq.init( depth );
	}
	int lca( int a , int b ){
		int x = in[a];
		int y = in[b];
		if( x > y ){
			swap( x , y );
		}
		int pos = rmq.query( x, y + 1 );
		return id[pos];
	}
};

LCA lca;
