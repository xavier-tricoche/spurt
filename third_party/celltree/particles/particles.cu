#include <glew.h>
#include <glut.h>

#include <dataset.hpp>
#include <celltree.hpp>
#include <celltree_builder.hpp>
#include <interpolator.hpp>
#include <cuda_interpolator.hpp>

#include "trackball.hpp"
#include <cstdio>
#include <cassert>
#include <cutil_math.h>
#include <cuda_gl_interop.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include "glm.h"
#include "fps_timer.hpp"

// -------------------------------------------------------------------------

interpolator*       intp_cpu;
cuda_interpolator*  intp_gpu;

trackball tb;
fps_timer fps;

const unsigned int N = 250;
const unsigned int M = 200;

const float T = 2.0;
const float H = 0.01;

float globalScale = 1.0;
float globalCenter[3] = { 0, 0, 0 };
bool  doAdvect    = true;
bool  advectGPU   = false;

GLuint contextModel = 0;
GLuint particleVBO  = 0;
cudaGraphicsResource* particleRes;

// -------------------------------------------------------------------------

struct particle
{
    float3 p;
    float  t;
    float4 c;
};

struct advect_cpu
{
    void operator()( particle& p ) const
    {
        p.t += H;

        if( p.t > 0 )
        {
            float3 v;

            if( !(*intp_cpu)( 0.0, &(p.p.x), &(v.x) ) )
                v = make_float3( 0, 0, 0 );

            float h = fmaxf( 0.0, fminf( p.t, H ) );

            p.p += h*v;
        }
    }
};

struct advect_gpu
{
    cuda_interpolator::interpolator intp;

    advect_gpu() : intp( intp_gpu->get_interpolator() )
    {
    }

    __device__
    void operator()( particle& p ) const
    {
        p.t += H;

        if( p.t > 0 )
        {
            float h = fmaxf( 0.0, fminf( p.t, H ) );

            p.p += h*intp( p.p );
        }
    }
};

// -------------------------------------------------------------------------

void particles_init()
{
    particle* p = new particle[N*M];

    unsigned int n = 0;

    // create initial particles; make sure that there is mostly
    // coherency between contiguous particles
    for( unsigned int iy=0; iy<M; ++iy )
    {
        for( unsigned int ix=0; ix<N; ++ix, ++n )
        {
            float rx = (ix -0.5 + drand48())/N;
            float ry = (iy -0.5 + drand48())/M;
            float rt = (iy -0.5 + drand48())/M;

            // ellipsoid
            float x =  -0.25 + 0.5*rx;
            float y =  -0.001 + 0.002*rt;
            float z =  0.3;
            float t = -T*rt;

            // tdelta
            // float x = -20;// 0.0  + 150.0*ry;
            // float y = -50.0 + 100.0*rx;
            // float z = -25.0; //-30.0 + 25.0*ry;
            // float t = -T*rt;

            // edelta
            // float x = 0.271009;
            // float y = 0.0865657 + 0.001*cos( 2.0*6.28318*rx );
            // float z = 0.0362513 + 0.001*sin( 2.0*6.28318*rx );
            // float t = -T*rt;

            // mx2all
            // float x = 0.0 + 5.5*cos( 2.0*6.28318*rx );
            // float y = 0.0 + 5.5*sin( 2.0*6.28318*rx );
            // float z = -10;
            // float t = -T*rt;

            // F6
            // float x = -190;
            // float y = -10 + rx*20;
            // float z = 94;
            // float t = -T*rt;

            // BMW
            // float x = -1900;
            // float y = -600+1200*rx;
            // float z = -200+450*rt;
            // float t = -T*rt;

            p[n].t = t;
            p[n].p = make_float3( x, y, z );
            p[n].c = make_float4( 1.0-0.9*rt, 0.5, 0.1+0.9*rt, 0.1 );
        }
    }

    // set up vertex buffer object
    if( particleVBO == 0 )
        glGenBuffers( 1, &particleVBO );

    glBindBuffer( GL_ARRAY_BUFFER, particleVBO );
    glBufferData( GL_ARRAY_BUFFER, sizeof(particle)*N*M, p, GL_STATIC_DRAW );
    glBindBuffer( GL_ARRAY_BUFFER, 0 );

    delete[] p;

    if( cudaSuccess != cudaGraphicsGLRegisterBuffer( &particleRes, particleVBO, cudaGraphicsMapFlagsNone ) )
        throw std::runtime_error( "unable to register VBO with CUDA" );

    glutReportErrors();
}

void particles_render()
{
    glPushAttrib( GL_ENABLE_BIT |
                  GL_COLOR_BUFFER_BIT |
                  GL_DEPTH_BUFFER_BIT );

    glDisable( GL_LIGHTING );
    glDepthMask( GL_FALSE );

    glBlendFunc( GL_SRC_ALPHA, GL_ONE );
    glEnable( GL_BLEND );
    glEnable( GL_POINT_SMOOTH );

    glPointSize( 1.1f );
    glColor4f( 1,0.5,0.1,0.1 );

    glBindBuffer( GL_ARRAY_BUFFER, particleVBO );

    glVertexPointer( 3, GL_FLOAT, sizeof(particle), (GLvoid*)0 );
    glColorPointer( 4, GL_FLOAT, sizeof(particle), (GLvoid*)16 );

    glEnableClientState( GL_VERTEX_ARRAY );
    glEnableClientState( GL_COLOR_ARRAY );

    glDrawArrays( GL_POINTS, 0, N*M );

    glDisableClientState( GL_VERTEX_ARRAY );
    glDisableClientState( GL_COLOR_ARRAY );
    glBindBuffer( GL_ARRAY_BUFFER, 0 );

    glPopAttrib();
}

void particles_advect_gpu()
{
    particle* p;
    size_t    s;

    if( cudaSuccess != cudaGraphicsMapResources( 1, &particleRes, 0 ) )
        throw std::runtime_error( "unable to map graphics resource" );

    if( cudaSuccess !=
        cudaGraphicsResourceGetMappedPointer( (void**)&p, &s, particleRes ) )
        throw std::runtime_error( "unable to obtain mapped resource pointer" );

    thrust::for_each( thrust::device_ptr<particle>(p),
                      thrust::device_ptr<particle>(p) + N*M,
                      advect_gpu() );

    if( cudaSuccess != cudaGraphicsUnmapResources( 1, &particleRes, 0 ) )
        throw std::runtime_error( "unable to unmap graphics resource" );
};

void particles_advect_cpu()
{
    glBindBuffer( GL_ARRAY_BUFFER, particleVBO );

    particle* p = (particle*)glMapBuffer( GL_ARRAY_BUFFER, GL_READ_WRITE );

    if( !p )
        throw std::runtime_error( "unable to map particle VBO" );

    std::for_each( p, p + N*M, advect_cpu() );

    glUnmapBuffer( GL_ARRAY_BUFFER );
    glBindBuffer( GL_ARRAY_BUFFER, 0 );
}

void initialize()
{
    particles_init();

    glClearColor( .1, .1, .1, 1.0 );

    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE );

    float specular[4] = { 0.5, 0.5, 0.5, 1.0 };
    glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, specular );
    glMaterialf( GL_FRONT_AND_BACK, GL_SHININESS, 50.0 );

    glEnable( GL_LIGHT0 );
}

void display()
{
    // count frame for fps count
    fps.frame();

    // overall setup
    glMatrixMode( GL_MODELVIEW );
    glLoadMatrixd( tb.matrix() );
    glScalef( globalScale, globalScale, globalScale );

    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    // draw context model
    glEnable( GL_LIGHTING );
    glDisable( GL_BLEND );
    glEnable( GL_DEPTH_TEST );
    glEnable( GL_LIGHTING );
    glEnable( GL_LIGHT0 );

    glCallList( contextModel );

    // draw particles
    particles_render();

    // output FPS to screen
    glDisable( GL_DEPTH_TEST );
    glDisable( GL_LIGHTING );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glMatrixMode( GL_PROJECTION );
    glPushMatrix();
    glLoadIdentity();

    glColor3f( 1, 1, 1 );
    glRasterPos3f( -0.97, -0.95, 0 );

    char buf[256];
    sprintf( buf, "%s FPS: %.0f", advectGPU ? "GPU" : "CPU", fps.rate() );

    for( char* c=buf; *c!=0; ++c )
        glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, *c );

    glPopMatrix();

    glutSwapBuffers();
    glutPostRedisplay();

    // do some advection
    if( doAdvect )
    {
        if( advectGPU )
            particles_advect_gpu();
        else
            particles_advect_cpu();
    }
}

void reshape( int w, int h )
{
    double aspect = (double)w / (double)h;

    const double fovy  = 45.0;
    const double zNear = 0.1;
    const double zFar  = 1000.0;

    double tan_fovy = tan( 0.5*fovy * M_PI/180.0 );
    double right    = tan_fovy * aspect * zNear;
    double left     = -right;
    double top      = tan_fovy * zNear;
    double bottom   = -top;

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glFrustum( left, right, bottom, top, zNear, zFar );

    glViewport( 0, 0, w, h );
}

void motion( int x, int y )
{
    tb.move( x/400.0 - 1.0, 1.0 - y/400.0 );
    glutPostRedisplay();
}

void mouse( int b, int state, int x, int y )
{
    int mod = glutGetModifiers();

    trackball::Mode mode = trackball::NONE;

    switch( b )
    {
    case GLUT_LEFT_BUTTON:
        mode = mod & GLUT_ACTIVE_SHIFT ? trackball::ROTATE_Z :
                                         trackball::ROTATE;
        break;
    case GLUT_RIGHT_BUTTON:
        mode = mod & GLUT_ACTIVE_SHIFT ? trackball::PUSH :
                                         trackball::ZOOM;
        break;
    case GLUT_MIDDLE_BUTTON:
        mode = trackball::PAN;
        break;
    }

    if( state == GLUT_DOWN )
        tb.begin( mode, x/400.0 - 1.0, 1.0 - y/400.0 );
    else
        tb.end();

    glutPostRedisplay();
}

void key( unsigned char key, int x, int y )
{
    int mod = glutGetModifiers();

    switch( key )
    {
    case 27:
        exit( 0 );
    case ' ':
        tb.home();
        break;
    case 'a':
        doAdvect = !doAdvect;
        break;
    case 'g':
        advectGPU = !advectGPU;
        break;
    case 'r':
        particles_init();
        break;
    default:
        std::cout << "unknown key " << (int)key << '\n';
        return;
    }

    glutPostRedisplay();
}

// ----------------------------------------------------

int main( int argc, char *argv[] )
{
    glutInit( &argc, argv );
    glutInitWindowSize( 960, 540 );

    glutInitDisplayMode( GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE );
    glutCreateWindow( argv[1] );

    glutReshapeFunc( reshape );
    glutDisplayFunc( display );
    glutMouseFunc( mouse );
    glutMotionFunc( motion );
    glutKeyboardFunc( key );

    glewInit();

    cudaSetDevice(0);
    cudaGLSetGLDevice(0);

    try
    {
        initialize();

        // load dataset and initialize interpolators
        dataset* ds = dataset::create( argv[1] );

        mesh* m     = ds->read_mesh();
        variable* v = ds->read_vector_variable( 0, "velocity" );

        celltree ct;
        celltree_builder builder;
        builder.build( ct, *m );

        intp_gpu = new cuda_interpolator( m, v, &ct );
        intp_cpu = new interpolator( m, v, ct );

        // load context model
        std::string prefix = argv[1];
        prefix = prefix.substr( 0, prefix.rfind( '.' ) );

        GLMmodel* obj = glmReadOBJ( (prefix+".obj").c_str() );

        glmBoundingSphere( obj, globalCenter, &globalScale );
        globalScale = 1.0 / globalScale;

        tb.home( 0, 0, 4, 0, 0, 0, 0, 1, 0 );

        glmFacetNormals( obj );
        contextModel = glmList( obj, GLM_FLAT );

        glmDelete( obj );

        // everything set up, let's go
        glutMainLoop();
    }
    catch( std::exception& e )
    {
        printf( "exception: %s\n", e.what() );
        exit( EXIT_FAILURE );
    }

    exit( EXIT_SUCCESS );
}
