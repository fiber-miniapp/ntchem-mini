#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <cublas_v2.h>

#define ERR_EQ(X,Y) do { if ((X) == (Y)) { \
    fprintf(stderr,"Error in %s at %s:%d\n",__func__,__FILE__,__LINE__); \
    exit(-1);}} while(0)

#define ERR_NE(X,Y) do { if ((X) != (Y)) { \
    fprintf(stderr,"Error in %s at %s:%d\n",__func__,__FILE__,__LINE__); \
    exit(-1);}} while(0)

#define CUDA_CALL(X)   ERR_NE(X,cudaSuccess)
#define CUBLAS_CALL(X) ERR_NE(X,CUBLAS_STATUS_SUCCESS)

// #define DPRINT(X, ...) fprintf( stderr, "[%s] " X, __func__, __VA_ARGS__ )
#define DPRINT(X, ...)

/*
 *
 */
#define MAX_STREAM  8
#define MAX_EVENT  (MAX_STREAM*2)

static cublasHandle_t  handle;
static cudaStream_t    stream[MAX_STREAM+1];
static cudaEvent_t     event [MAX_EVENT+1];
static cudaEvent_t     event2[MAX_EVENT+1];

static int  init_flag = 0;

extern "C"
void cublas_init_()
{
    if ( init_flag ) return;

    CUBLAS_CALL( cublasCreate( &handle ) );
    for ( int i = 0; i <= MAX_STREAM; i++ ) {
	CUDA_CALL( cudaStreamCreate( &stream[i] ) );
    }
    for ( int i = 0; i <= MAX_EVENT; i++ ) {
	CUDA_CALL( cudaEventCreate( &event[i] ) );
	CUDA_CALL( cudaEventCreate( &event2[i] ) );
    }

    init_flag = 1;
}

extern "C"
void cublas_fin_()
{
    if ( ! init_flag ) return;

    for ( int i = 0; i <= MAX_EVENT; i++ ) {
	CUDA_CALL( cudaEventDestroy( event[i] ) );
	CUDA_CALL( cudaEventDestroy( event2[i] ) );
    }
    for ( int i = 0; i <= MAX_STREAM; i++ ) {
	CUDA_CALL( cudaStreamDestroy( stream[i] ) );
    }
    CUBLAS_CALL( cublasDestroy( handle ) );

    init_flag = 0;
}

extern "C"
void cublas_alloc_( void **devptr, const int *i, const int *j )
{
    size_t  size = sizeof(double) * i[0] * j[0];
    DPRINT( "size:%ld (%ld x %ld x %ld)\n", size, sizeof(double), i[0], j[0] );
    void  *buf;
    CUDA_CALL( cudaMalloc( &buf, size ) );
    *devptr = buf;
    DPRINT( "devptr:%p, size:%ld (%ld x %ld x %ld)\n", *devptr, size, sizeof(double), i[0], j[0] );
}

extern "C"
void cublas_free_( void **devptr )
{
    DPRINT( "devptr:%p\n", *devptr );
    void  *buf;
    buf = *devptr;
    CUDA_CALL( cudaFree( buf ) );
}

extern "C"
void cublas_set_matrix_( const int *i, const int *j, const void *A, const int *lda,
		      void **devptr_B, const int *ldb )
{
    DPRINT( "A:%p, devptr_B:%p\n", A, *devptr_B );
    double  *dev_B;
    dev_B = (double*) *devptr_B;
    CUBLAS_CALL( cublasSetMatrix( *i, *j, sizeof(double), A, *lda, dev_B, *ldb ) );
}

extern "C"
void cublas_set_matrix_async_( const int *i, const int *j, const void *A, const int *lda,
			       void **devptr_B, const int *ldb, const int *id_st )
{
    DPRINT( "A:%p, devptr_B:%p, st:%d\n", A, *devptr_B, *id_st );
    assert( *id_st >= 0 && *id_st <= MAX_STREAM );
    double  *dev_B;
    dev_B = (double*) *devptr_B;
    CUBLAS_CALL( cublasSetMatrixAsync( *i, *j, sizeof(double), A, *lda, dev_B, *ldb, stream[*id_st] ) );
}

extern "C"
void cublas_get_matrix_( const int *i, const int *j, const void **devptr_A, const int *lda, 
		      void *B, const int *ldb )
{
    DPRINT( "devptr_A:%p, B:%p\n", *devptr_A, B );
    double  *dev_A;
    dev_A = (double*) *devptr_A;
    CUBLAS_CALL( cublasGetMatrix( *i, *j, sizeof(double), dev_A, *lda, B, *ldb ) );
}

extern "C"
void cublas_get_matrix_async_( const int *i, const int *j, const void **devptr_A, const int *lda, 
			       void *B, const int *ldb, const int *id_st )
{
    DPRINT( "devptr_A:%p, B:%p, st:%d\n", *devptr_A, B, *id_st );
    assert( *id_st >= 0 && *id_st <= MAX_STREAM );
    double  *dev_A;
    dev_A = (double*) *devptr_A;
    CUBLAS_CALL( cublasGetMatrixAsync( *i, *j, sizeof(double), dev_A, *lda, B, *ldb, stream[*id_st] ) );
}

extern "C"
void cublas_dgemm_( const char *ta, const char *tb,
		 const int *m, const int *n, const int *k,
		 const double *alpha,
		 const void **devptr_A, const int *lda,
		 const void **devptr_B, const int *ldb,
		 const double *beta,
		 void **devptr_C, const int *ldc )
{
    DPRINT( "m:%d, n:%d, k:%d\n", *m, *n, *k );
    cublasOperation_t  transa = CUBLAS_OP_N;
    cublasOperation_t  transb = CUBLAS_OP_N;
    if ( toupper(ta[0]) == 'T' ) transa = CUBLAS_OP_T;
    if ( toupper(tb[0]) == 'T' ) transb = CUBLAS_OP_T;
    const double  *dev_A;
    const double  *dev_B;
    double  *dev_C;
    dev_A = (double*) *devptr_A;
    dev_B = (double*) *devptr_B;
    dev_C = (double*) *devptr_C;
    CUBLAS_CALL ( cublasDgemm( handle, transa, transb, *m, *n, *k, alpha,
			       dev_A, *lda,
			       dev_B, *ldb, beta,
			       dev_C, *ldc ) );
}

extern "C"
void cublas_dgemm_async_( const char *ta, const char *tb,
			  const int *m, const int *n, const int *k,
			  const double *alpha,
			  const void **devptr_A, const int *lda,
			  const void **devptr_B, const int *ldb,
			  const double *beta,
			  void **devptr_C, const int *ldc,
			  const int *id_st )
{
    DPRINT( "m:%d, n:%d, k:%d, st:%d\n", *m, *n, *k, *id_st );
    assert( *id_st >= 0 && *id_st <= MAX_STREAM );
    cublasOperation_t  transa = CUBLAS_OP_N;
    cublasOperation_t  transb = CUBLAS_OP_N;
    if ( toupper(ta[0]) == 'T' ) transa = CUBLAS_OP_T;
    if ( toupper(tb[0]) == 'T' ) transb = CUBLAS_OP_T;
    const double  *dev_A;
    const double  *dev_B;
    double  *dev_C;
    dev_A = (double*) *devptr_A;
    dev_B = (double*) *devptr_B;
    dev_C = (double*) *devptr_C;
    CUBLAS_CALL ( cublasSetStream( handle, stream[*id_st] ) );
    CUBLAS_CALL ( cublasDgemm( handle, transa, transb, *m, *n, *k, alpha,
			       dev_A, *lda,
			       dev_B, *ldb, beta,
			       dev_C, *ldc ) );
    CUBLAS_CALL ( cublasSetStream( handle, NULL ) );
}

extern "C"
void cublas_st_sync_( const int *id_st )
{
    DPRINT( "st:%d\n", *id_st );
    assert( *id_st >= 0 && *id_st <= MAX_STREAM );
    CUDA_CALL( cudaStreamSynchronize( stream[*id_st] ) );
}

extern "C"
void cublas_ev_rec_( const int *id_ev, const int *id_st )
{
    DPRINT( "ev:%d, st:%d\n", *id_ev, *id_st );
    assert( *id_st >= 0 && *id_st <= MAX_STREAM );
    assert( *id_ev >= 0 && *id_ev <= MAX_EVENT );
    CUDA_CALL( cudaEventRecord( event[*id_ev], stream[*id_st] ) );
}

extern "C"
void cublas_ev_wait_( const int *id_ev, const int *id_st )
{
    DPRINT( "ev:%d, st:%d\n", *id_ev, *id_st );
    assert( *id_st >= 0 && *id_st <= MAX_STREAM );
    assert( *id_ev >= 0 && *id_ev <= MAX_EVENT );
    CUDA_CALL( cudaStreamWaitEvent( stream[*id_st], event[*id_ev], 0 ) );
}

extern "C"
void cublas_ev2_rec_( const int *id_ev, const int *id_st )
{
    DPRINT( "ev:%d, st:%d\n", *id_ev, *id_st );
    assert( *id_st >= 0 && *id_st <= MAX_STREAM );
    assert( *id_ev >= 0 && *id_ev <= MAX_EVENT );
    CUDA_CALL( cudaEventRecord( event2[*id_ev], stream[*id_st] ) );
}

extern "C"
void cublas_ev2_wait_( const int *id_ev, const int *id_st )
{
    DPRINT( "ev:%d, st:%d\n", *id_ev, *id_st );
    assert( *id_st >= 0 && *id_st <= MAX_STREAM );
    assert( *id_ev >= 0 && *id_ev <= MAX_EVENT );
    CUDA_CALL( cudaStreamWaitEvent( stream[*id_st], event2[*id_ev], 0 ) );
}
