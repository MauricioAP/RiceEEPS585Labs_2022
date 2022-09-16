/*
 *
 * Author: Mauricio Araya Polo
 * Date: 04/2015 - present
 *
 */

#include "props2D.h"

void absmult( float *a, float *b, float *c, int *prms )
{
    uint iz, ix, column, indxe;
    int stencil = prms[0];
    int nz = prms[1];
    int nx = prms[2];
    int col = prms[9];

    #pragma omp parallel for private(iz, ix, column, indxe), schedule(runtime)
    for( ix = stencil; ix < nx + stencil; ix++ )
    {
        column = ix*col;
        for( iz = stencil + 1; iz < nz + stencil; iz++ )
        {
            indxe = iz + column;
            a[ indxe ] *= c[ indxe ];
            b[ indxe ] *= c[ indxe ];

	}
    }
}

void mmult( float *a, float *b, float *c, int *prms, float *d, float *e )
{
    uint iz, ix, column, indxe;
    int stencil = prms[0];
    int nz = prms[1];
    int nx = prms[2];
    int col = prms[9];

    #pragma omp parallel for private(iz, ix, column, indxe), schedule(runtime)
    for( ix = stencil; ix < nx + stencil; ix++ )
    {
        column = ix*col;
        for( iz = stencil + 1; iz < nz + stencil; iz++ )
        {
            indxe = iz + column;
            c[ indxe ] += a[ indxe ] * b[ indxe ];
            d[ indxe ] += b[ indxe ] * b[ indxe ];
            e[ indxe ] += a[ indxe ] * a[ indxe ];

	}
    }
}

void fwd_step2D( float *p, float *po, float *attr, int *prms, float *cfs )
{
    float tmp2, tmp3[4000];
    uint indx, column, indxe, iz, ix;

    TRANSLATE2D()

    #pragma omp parallel for private(iz, ix, indx, indxe, tmp2, tmp3, column), schedule(runtime)
    for( ix = stencil; ix < nx + stencil; ix++ )
    {
        column = ix*col;
        // manual pre-fetching, non-contiguous memory access
        indxe = stencil + column;
        tmp3[0]  = c0 * po[ indxe ] +
                   cx1 * ( po[ indxe + col  ] + po[ indxe - col  ] ) +
                   cx2 * ( po[ indxe + col2 ] + po[ indxe - col2 ] ) +
                   cx3 * ( po[ indxe + col3 ] + po[ indxe - col3 ] ) +
                   cx4 * ( po[ indxe + col4 ] + po[ indxe - col4 ] ) ;
        #pragma omp simd
        for( iz = stencil + 1; iz < nz + stencil; iz++ )
        {
            indxe = iz + column;

            tmp3[iz-stencil] = c0 * po[ indxe ] +
                               cx1 * ( po[ indxe + col  ] + po[ indxe - col  ] ) +
                               cx2 * ( po[ indxe + col2 ] + po[ indxe - col2 ] ) +
                               cx3 * ( po[ indxe + col3 ] + po[ indxe - col3 ] ) +
                               cx4 * ( po[ indxe + col4 ] + po[ indxe - col4 ] ) ;
        }
        #pragma omp simd
        for( iz = stencil; iz < nz + stencil; iz++ )
        {
            indxe = iz + column;
            indx = (iz-stencil) + (ix-stencil) * nz;
            tmp2 = attr[ indx ] * dt;
            tmp2 *= tmp2;

            tmp3[iz-stencil] = tmp2*( tmp3[iz-stencil] +
                               cz1 * ( po[ indxe + 1 ] + po[ indxe - 1 ] ) +
                               cz2 * ( po[ indxe + 2 ] + po[ indxe - 2 ] ) +
                               cz3 * ( po[ indxe + 3 ] + po[ indxe - 3 ] ) +
                               cz4 * ( po[ indxe + 4 ] + po[ indxe - 4 ] )
                               );

            p[ indxe ] = po[ indxe ] + po[ indxe ] + tmp3[iz-stencil] - p[ indxe ];

        }
    }
}
