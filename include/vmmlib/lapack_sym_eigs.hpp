#ifndef __VMML__VMMLIB_LAPACK_SYM_EIGS__HPP__
#define __VMML__VMMLIB_LAPACK_SYM_EIGS__HPP__

#include <vmmlib/matrix.hpp>
#include <vmmlib/vector.hpp>
#include <vmmlib/exception.hpp>

#include <vmmlib/lapack_types.hpp>
#include <vmmlib/lapack_includes.hpp>

#include <string>

/** 
 *
 *   a wrapper for lapack's DSYEVX routine. 
 *   DSYEVX computes selected eigenvalues and, optionally, eigenvectors
 *   of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
 *   selected by specifying either a range of values or a range of indices
 *   for the desired eigenvalues.
 *   
 *   returns a boolean to indicate success of the operation. 
 *   if the return value is false, you can get the parameters using
 *   get_params(). 
 *   error states:
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit.
 *          < 0:  if INFO = -i, the i-th argument had an illegal value.
 *          > 0:  if DBDSQR did not converge, INFO specifies how many
 *                superdiagonals of an intermediate bidiagonal form B
 *                did not converge to zero. See the description of WORK
 *                above for details.
 *
 *   more information in: http://www.netlib.org/lapack/double/dsyevx.f (see also:
 *   http://www.netlib.org/lapack/double/dsyev.f , but    needs more space)
 **
 */

namespace vmml
{
	
namespace lapack
{
	
	// XYYZZZ 
	// X    = data type: S - float, D - double
	// YY   = matrix type, GE - general, TR - triangular
	// ZZZ  = function name
	
	
	template< typename float_t >
	struct eigs_params
	{
		char            jobz;
		char            range;
		char            uplo;
		lapack_int      n;
		float_t*        a;
		lapack_int      lda; //leading dimension of input array
		float_t*        vl;
		float_t*        vu;
		lapack_int      il;
		lapack_int      iu;
		float_t         abstol;
		lapack_int      m; //number of found eigenvalues
		float_t*        w; //first m eigenvalues
		float_t*        z; //first m eigenvectors
		lapack_int      ldz; //leading dimension of z
		float_t*        work;
		lapack_int      lwork;
		lapack_int*     iwork;
		lapack_int*     ifail;
		lapack_int      info;
		
		friend std::ostream& operator << ( std::ostream& os, 
										  const eigs_params< float_t >& p )
		{
			os 
			<< " (1)\tjobz "     << p.jobz << std::endl
			<< " (2)\trange "    << p.range << std::endl
			<< " (3)\tuplo "     << p.uplo << std::endl
			<< " (4)\tn "        << p.n << std::endl
			<< " (5)\ta "        << *p.a << std::endl
			<< " (6)\tlda "      << p.lda << std::endl
			<< " (7)\tvl "       << p.vl << std::endl 
			<< " (8)\tvu "       << p.vu << std::endl
			<< " (9)\til "       << p.il << std::endl
			<< " (10)\tiu "       << p.iu << std::endl
			<< " (11)\tabstol "   << p.abstol << std::endl
			<< " (12)\tm "        << p.m << std::endl
			<< " (13)\tw "        << p.w << std::endl
			<< " (14)\tz "        << p.z << std::endl
			<< " (15)\tldz "      << p.ldz << std::endl
			<< " (16)\twork "     << *p.work << std::endl
			<< " (17)\tlwork "    << p.lwork << std::endl
			<< " (18)\tiwork "    << *p.iwork << std::endl
			<< " (19)\tifail "    << *p.ifail << std::endl
			<< " (20)\tinfo "     << p.info 
			<< std::endl;
			return os;
		}
		
	};
	
	
#if 0
	/* Subroutine */ 

	int dsyevx_( char *jobz, char *range, char *uplo, integer *n, doublereal *a, integer *lda, doublereal *vl, 
				doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *w, doublereal *z, 
				integer *ldz, doublereal *work, integer* lwork, integer* iwork, integer* ifail, integer* info );
	
#endif
	
	
	template< typename float_t >
	inline void
	sym_eigs_call( eigs_params< float_t >& p )
	{
		VMMLIB_ERROR( "not implemented for this type.", VMMLIB_HERE );
	}
	
	
	template<>
	inline void
	sym_eigs_call( eigs_params< float >& p )
	{
		//std::cout << "calling lapack sym x eigs (single precision) " << std::endl;
		ssyevx_( 
				&p.jobz,
				&p.range,
				&p.uplo,
				&p.n,
				p.a,
				&p.lda,
				p.vl,
				p.vu,
				&p.il,
				&p.iu,
				&p.abstol,
				&p.m,
				p.w,
				p.z,
				&p.ldz,
				p.work,
				&p.lwork,
				p.iwork,
				p.ifail,
				&p.info
				);
		
	}
	
	
	template<>
	inline void
	sym_eigs_call( eigs_params< double >& p )
	{
		//std::cout << "calling lapack sym x eigs (double precision) " << std::endl;
		dsyevx_( 
				&p.jobz,
				&p.range,
				&p.uplo,
				&p.n,
				p.a,
				&p.lda,
				p.vl,
				p.vu,
				&p.il,
				&p.iu,
				&p.abstol,
				&p.m,
				p.w,
				p.z,
				&p.ldz,
				p.work,
				&p.lwork,
				p.iwork,
				p.ifail,
				&p.info
				);
	}
	
} // namespace lapack



template< size_t N, typename float_t >
struct lapack_sym_eigs
{
	lapack_sym_eigs();
	~lapack_sym_eigs();
	
	// version of reduced sym. eigenvalue decomposition, 
	// computes only the x largest magn. eigenvalues and their corresponding eigenvectors 
	
	//computes all eigenvalues and eigenvectors for matrix A
	bool compute_all(
				 const matrix< N, N, float_t >& A,
				 matrix< N, N, float_t >& eigvectors,
				 vector< N, float_t >& eigvalues
				 );
		
	inline bool test_success( lapack::lapack_int info );
	
	lapack::eigs_params< float_t > p;
	
	const lapack::eigs_params< float_t >& get_params(){ return p; };
	
}; // struct lapack_sym_eigs


template< size_t N, typename float_t >
lapack_sym_eigs< N, float_t >::lapack_sym_eigs()
{
	p.jobz      = 'V'; // Compute eigenvalues and eigenvectors.
	p.range     = 'A'; // all eigenvalues will be found.
	p.uplo      = 'U'; // Upper triangle of A is stored; or Lower triangle of A is stored.
	p.n         = N;
	p.a         = 0; //If UPLO = 'U', the leading N-by-N upper triangular part of A contains the upper triangular part of the matrix A.
	p.lda       = N;
	p.vl        = 0; //Not referenced if RANGE = 'A' or 'I'.
	p.vu        = 0; //Not referenced if RANGE = 'A' or 'I'.	 
	p.il        = 0; //Not referenced if RANGE = 'A' or 'V'.
	p.iu        = 0; //Not referenced if RANGE = 'A' or 'V'.	
	p.abstol    = 0.0001; //lie in an interval [a,b] of width less than or equal to ABSTOL + EPS *   max( |a|,|b| )
	p.m         = N; //The total number of eigenvalues found.  0 <= M <= N.
	p.w         = 0; //first m eigenvalues
	p.z         = 0; //first m eigenvectors
	p.ldz       = N; // The leading dimension of the array Z.  LDZ >= 1, and if JOBZ = 'V', LDZ >= max(1,N).
	p.work      = new float_t;
	//FIXME: check if correct datatype
	p.iwork     = new __CLPK_integer[5*N]; //[5*N]; // INTEGER array, dimension (5*N)
	p.ifail     = new __CLPK_integer[N]; //[N];
	p.lwork     = -1; //8N
	
	// workspace query
	lapack::sym_eigs_call( p );
	
	p.lwork = static_cast< lapack::lapack_int >( p.work[0] );
	delete p.work;
	
	p.work = new float_t[ p.lwork ];
	
}



template< size_t N, typename float_t >
lapack_sym_eigs< N, float_t >::~lapack_sym_eigs()
{
	delete[] p.work; 
	delete[] p.iwork;
	delete[] p.ifail;
}
	
template< size_t N, typename float_t >
bool
lapack_sym_eigs< N, float_t >::compute_all(
									   const matrix< N, N, float_t >& A,
									   matrix< N, N, float_t >& eigvectors,
									   vector< N, float_t >& eigvalues
									   )
{
	// lapack destroys the contents of the input matrix
	matrix< N, N, float_t > AA( A );

	p.range     = 'A'; // all eigenvalues will be found.
	p.a         = AA.array;
	p.ldz       = N;
    p.w         = eigvalues.array;
    p.z         = eigvectors.array;
	
	//debug std::cout << p << std::endl;
	
	lapack::sym_eigs_call< float_t >( p );
	
	return p.info == 0;
}	

	
	
} // namespace vmml

#endif	
