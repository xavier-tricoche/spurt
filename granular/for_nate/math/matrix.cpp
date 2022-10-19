#include "math/matrix.hpp"
#include "math/math.hpp"

bool xavier::display_eigen_stuff;

// ---------------------------------------------------------------------------

void xavier::mat3::eigensystem( std::vector< double >& evals, 
	std::vector< nvis::vec3 >& evecs ) const
{
	// coefficients of characteristic polynomial
	double coef[4] = { -1, __mat[0]+__mat[4]+__mat[8],
		__mat[6]*__mat[2]+__mat[7]*__mat[5]+
		__mat[3]*__mat[1]-__mat[0]*__mat[4]-
		__mat[0]*__mat[8]-__mat[4]*__mat[8],
		__mat[0]*__mat[4]*__mat[8]+
		__mat[1]*__mat[5]*__mat[6]+
		__mat[2]*__mat[3]*__mat[7]-
		__mat[6]*__mat[4]*__mat[2]-
		__mat[0]*__mat[5]*__mat[7]-
		__mat[8]*__mat[3]*__mat[1] };

	// solve cubic equation
	std::complex< double > __evals[3]; 
	int n = cubic_equation( coef[0], coef[1], coef[2], coef[3], __evals );

	// find corresponding eigenvectors
	evecs.clear();
	evals.clear();
	for ( unsigned int i=0 ; i<n ; i++ )
	{
		if ( __evals[i].imag() ) continue;

		double lambda = __evals[i].real();
		nvis::vec3 ev = eigenvector( *this, lambda );
		if ( display_eigen_stuff && nvis::norm( ev )==0 )
		{
			std::cout << "unable to find an eigenvector for eigenvalue " 
				<< lambda << std::endl;
			continue;
		}

		evals.push_back( __evals[i].real() );
		evecs.push_back( ev );

		if ( display_eigen_stuff )
		{
			// checking solution
			std::cout << "checking eigenvalue #" << i << ": " << lambda << std::endl
				<< "on characteristic polynomial: " 
				<< -lambda*lambda*lambda + coef[1]*lambda*lambda + coef[2]*lambda + coef[3]
				<< std::endl
				<< "on matrix: " << std::flush;
			mat3 __M( *this );
			__M(0,0) -= lambda;
			__M(1,1) -= lambda;
			__M(2,2) -= lambda;
			std::cout << det(__M) << std::endl;

			std::cout << "checking eigenvector #" << i << ": " << ev << std::endl;
			nvis::vec3 Me = __M*ev;
			std::cout << "(M-lambda I_3)*ev = " << Me << std::endl;
			Me = (*this)*ev;
			std::cout << "< M*ev, ev >/lambda = "
				<< nvis::inner( Me, ev )/lambda << std::endl;
			std::cout << "(M*ev) x ev = " << nvis::cross( Me, ev ) 
				<< std::endl;
		}
	}
}

// ---------------------------------------------------------------------------

void xavier::check_eigen( const mat3& M )
{
	if ( true || !display_eigen_stuff ) return;

	double m[9] = { M(0,0), M(0,1), M(0,2),
		M(1,0), M(1,1), M(1,2),
		M(2,0), M(2,1), M(2,2) };

	std::cout << "checking " << M << std::endl;

	double val[3];
	double vec[9];
	ell_3m_eigensolve_d( val, vec, m, 0 );

	std::vector< double > vals(3);
	std::vector< nvis::vec3 > vecs(3);
	for ( unsigned int i=0 ; i<3 ; i++ )
	{
		vals[i] = val[i];
		vecs[i] = nvis::vec3( vec[3*i], vec[3*i+1], vec[3*i+2] );
	}

	if ( vals[0] == vals[1] || vals[1] == vals[2] )
	{
		std::cout << "redundant eigenvalues: " 
			<< vals[0] << " >= " << vals[1] << " >= " << vals[2]
			<< std::endl;
	}

	for ( unsigned int i=0 ; i<3 ; i++ )
	{
		std::cout << "checking eigenvalue #" << i << ": " << vals[i] << std::endl;
		mat3 tmp( M );
		tmp(0,0) -= vals[i];
		tmp(1,1) -= vals[i];
		tmp(2,2) -= vals[i];

		std::cout << "det( M-lambda I_3 ) = " << det( tmp ) 
			<< " (det M = " << det(M) << ")" << std::endl;

		std::cout << "checking eigenvector #" << i << ": " << vecs[i] << std::endl;
		nvis::vec3 zero = tmp*vecs[i];
		std::cout << "( M-lambda I3 )*ev = " << zero << std::endl;

		std::cout << "M-lambda I_3 = " << tmp << std::endl;
		nvis::vec3 row0( tmp(0,0), tmp(0,1), tmp(0,2) );
		nvis::vec3 row1( tmp(1,0), tmp(1,1), tmp(1,2) );
		nvis::vec3 row2( tmp(2,0), tmp(2,1), tmp(2,2) );

		std::cout << "checking colinearity" << std::endl;
		double dot01 = fabs( nvis::inner( row0, row1 ) )/nvis::norm( row0 )/nvis::norm( row1 );
		double dot02 = fabs( nvis::inner( row0, row2 ) )/nvis::norm( row0 )/nvis::norm( row2 );
		double dot12 = fabs( nvis::inner( row1, row2 ) )/nvis::norm( row1 )/nvis::norm( row2 );
		std::cout << "0-1: " << dot01 << std::endl; 
		std::cout << "0-2: " << dot02 << std::endl;
		std::cout << "1-2: " << dot12 << std::endl;

		double s, t, delta;
		if ( dot01<dot02 && dot01<dot12 )
		{
			std::cout << "solving for 0-1" << std::endl;
			delta = row0[0]*row1[1]-row0[1]*row1[0];
			s = ( -row0[2]*row1[1]+row1[2]*row0[1] )/delta;
			t = ( -row0[0]*row1[2]+row1[0]*row0[2] )/delta;
		}
		else if ( dot02<dot01 && dot02<dot12 )
		{
			std::cout << "solving for 0-2" << std::endl;
			delta = row0[0]*row2[1]-row0[1]*row2[0];
			s = ( -row0[2]*row2[1]+row2[2]*row0[1] )/delta;
			t = ( -row0[0]*row2[2]+row2[0]*row0[2] )/delta;
		}
		else
		{
			std::cout << "solving for 1-2" << std::endl;
			delta = row1[0]*row2[1]-row1[1]*row2[0];
			s = ( -row1[2]*row2[1]+row2[2]*row1[1] )/delta;
			t = ( -row1[0]*row2[2]+row2[0]*row1[2] )/delta;
		}

		std::cout << "local coordinates are: (" << s << ", " << t << ")" << std::endl;
		nvis::vec3 q( s, t, 1 );
		q = tmp*q;
		std::cout << "( M-lambda I3 )*(s, t, 1)^T = " << q << std::endl;
	}
}

// ---------------------------------------------------------------------------

nvis::vec3 xavier::eigenvector( const mat3& M, const double lambda )
{
	static unsigned int indices[3][2]={ {1,2}, {0,2}, {0,1} };
	static const double epsilon = 1.0e-6;
	unsigned int *col, *row; 
	double det[3];
	unsigned int maxi(0);
	double maxdet;

		// M-lambda I_3
	mat3 matrix(M);
	double _norm = norm(matrix);
	for ( unsigned int j=0 ; j<3 ; ++j )
		matrix(j,j) -= lambda;

		// loop over scalar components of the eigenvector to 
		// find one that is not zero
	for ( unsigned int j=0 ; j<3 ; ++j )
	{
		col = indices[j];

			// assume evec[j]=1 i.e. is non zero
			// compute determinant of associated subsystems
		for ( unsigned int k=0 ; k<3 ; ++k )
		{
			row = indices[k];
			det[k] = matrix(row[0],col[0])*matrix(row[1],col[1])-
				matrix(row[1],col[0])*matrix(row[0],col[1]);
		}

			// get maximum
		maxdet = 0;
		for ( unsigned int k=0 ; k<3 ; ++k )
		{
			if ( fabs(det[k]) > maxdet )
			{
				maxi = k;
				maxdet = fabs(det[k]);
			}
		}
		if ( maxdet < epsilon*_norm )
			continue;

			// can solve the system by setting maxi-th coordinate equal to 1
		row = indices[maxi];
		det[maxi] = 1./det[maxi];

		nvis::vec3 evec;
		evec[col[0]] = det[maxi]*
			( -matrix( row[0], j )*matrix( row[1], col[1] )+
			matrix( row[1], j )*matrix( row[0], col[1] ) );
		evec[col[1]] = det[maxi]*
			( -matrix( row[0], col[0] )*matrix( row[1], j )+
			matrix( row[1], col[0] )*matrix( row[0], j ) );
		evec[j] = 1.;
		evec *= 1/nvis::norm( evec );

		return evec;
	}

	return nvis::vec3( 0, 0, 0 );
}