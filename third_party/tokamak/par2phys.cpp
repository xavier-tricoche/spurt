#include <cmath>
#include <iostream>

#include <GLUT/glut.h>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <boost/program_options.hpp>

#include "tokamak_nimrod_parametric.hpp"

#include <teem/nrrd.h>

using namespace nvis;

bbox2 bounds;
GLuint list, tex;

unsigned int rx, ry;
bool saveonly;
std::string nimrodfile;
std::string input;
std::string output;

float scale = 1.0;
float dx = 0.0, dy = 0.0;

int button, xstart, ystart;

void init()
{
    std::vector<vec2> parm, phys;
    unsigned int n1, n2;

    {
        tokamak_nimrod_parametric tnp( nimrodfile.c_str(), "15000" );
        tnp.get_slice( parm, phys, n1, n2 );
    }
    
    list = glGenLists(1);
    glNewList( list, GL_COMPILE );

    glBegin(GL_QUADS);

    unsigned int off[4] = { 0, 1, n1+1, n1 };

    for( unsigned int i=0; i<n2-1; ++i )
        for( unsigned int j=0; j<n1-1; ++j )
        {
            unsigned int base = n1*i+j;

            for( unsigned int k=0; k<4; ++k )
            {
                glTexCoord2f( parm[base+off[k]][0], parm[base+off[k]][1] );
                glVertex2f( phys[base+off[k]][0], phys[base+off[k]][1] );
            }
        }

    glEnd();

    glEndList();

    bounds = bbox2( phys.begin(), phys.end() );

    // ---
    
    Nrrd* nrrd = nrrdNew();
    
    if( nrrdLoad( nrrd, input.c_str(), 0 ) )
    {
        char *err = biffGetDone( NRRD );
        
        std::cerr << "trouble reading NRRD/PNG:" << err << '\n';
        exit( -1 );
    }
    
    Nrrd* cropped = nrrdNew();
    
    size_t min[] = { 0, 1, 1 }, max[] = { 2, nrrd->axis[1].size-2, nrrd->axis[2].size-2 };
    nrrdCrop( cropped, nrrd, min, max );
    nrrdNuke( nrrd );

    glGenTextures( 1, &tex );
    glBindTexture( GL_TEXTURE_2D, tex );
    
    glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );

    GLenum datatype;

    switch( cropped->type )
    {
    case nrrdTypeUChar:  datatype = GL_UNSIGNED_BYTE; break;
    case nrrdTypeUShort: datatype = GL_UNSIGNED_SHORT; break;
    default:
        std::cerr << "wrong format NRRD/PNG, need 3-dimensional uchar data\n";
        exit( -1 );
    }    

    gluBuild2DMipmaps( GL_TEXTURE_2D, GL_RGB, cropped->axis[1].size, cropped->axis[2].size,
                  GL_RGB, datatype, cropped->data );
                  
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 2 );
    
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
    
    glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );

    nrrdNuke( cropped );

    glutReportErrors();
}

void reshape( int w, int h )
{
    rx = w;
    ry = h;
}


void save()
{
    Nrrd* out = nrrdNew();
    nrrdAlloc_va( out, nrrdTypeUChar, 3, 4, rx, ry );

    glReadPixels( 0, 0, rx, ry, GL_RGBA, GL_UNSIGNED_BYTE, out->data );
    
    nrrdSave( output.c_str(), out, 0 );
}

void display(void)
{    
    glViewport( 0, 0, rx, ry );
    
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glOrtho( bounds.min()[0], bounds.max()[0], bounds.min()[1], bounds.max()[1], 0, 1 );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    glTranslatef( bounds.center()[0], bounds.center()[1], 0 );
    glScalef( scale, scale, 0.0 );
    glTranslatef( -bounds.center()[0], -bounds.center()[1], 0 );
    glTranslatef( dx, dy, 0 );
    
    glClearColor( 0.5, 0.5, 0.5, 0.0 );
    glClear( GL_COLOR_BUFFER_BIT );
    
    glPointSize(2);
    
    glEnable( GL_TEXTURE_2D );
    glCallList( list );
    
    glFinish();
    
    glutReportErrors();
}

void mouse( int b, int state, int x, int y )
{
    button = b;
    xstart = x;
    ystart = y;
}

void motion( int x, int y )
{
    if( button == GLUT_RIGHT_BUTTON )
        scale *= pow( 1.01, ystart-y );
    else if( button == GLUT_LEFT_BUTTON )
    {
        dx += (x-xstart)/float(rx) / scale;
        dy -= (y-ystart)/float(ry) / scale;
    }
    xstart = x;
    ystart = y;
    
    glutPostRedisplay();
}

void key( unsigned char key, int x, int y )
{
    switch( key )
    {
    case 27:
        exit( 0 );
    case ' ':
        scale = 1.0;
        dx = dy = 0.0;
        break;
    case 's':
        display();
        save();
        break;
    }

    glutPostRedisplay();
}

int main( int argc, char **argv )
{
    using namespace boost::program_options;
    
    options_description args( "Options" );
    
    args.add_options()
        ("help",    "produce help message")
        ("save-only,s", value<bool>( &saveonly )->default_value( false )->zero_tokens(), "save to file and quit immediately")
        ("rx", value<unsigned int>( &rx )->default_value( 768 ), "render resolution x" )
        ("ry", value<unsigned int>( &ry )->default_value( 768 ), "render resolution y" );
    
    args.add_options()
        ("nimrodfile",   value<std::string>( &nimrodfile ), "NIMROD HDF5 file" )
        ("input", value<std::string>( &input ), "input nrrd/png" )
        ("output", value<std::string>( &output ), "output nrrd/png" );
    
    positional_options_description pargs;
    pargs.add( "nimrodfile", 1 );
    pargs.add( "input",  1 );
    pargs.add( "output",  1 );

    std::ostringstream usage;
     
    usage << "Usage: " << argv[0] << " [options] <map description> <output nrrd>\n"
          << args << '\n';
    
    variables_map vm;
        
    try
    {
        store( command_line_parser( argc, argv ).options(args).positional(pargs).run(), vm );
        notify( vm );    
    }
    catch( std::exception& e )
    {
        std::cerr << "error parsing command line\n" << e.what() << "\n"
                  << usage.str();
        return -1;
    }
          
    if( vm.count("help") ) 
    {
        std::cerr << usage.str() << "\n";
        return -1;
    }

    if( !vm.count("nimrodfile") )
    {   
        std::cerr << "required nimrod file not given\n\n" << usage.str() << '\n';
        return -1;
    }
    
    if( !vm.count("input") )
    {   
        std::cerr << "input nrrd/png not given\n\n" << usage.str() << '\n';
        return -1;
    }
    
    if( !vm.count("output") )
    {   
        std::cerr << "output nrrd/png not given\n\n" << usage.str() << '\n';
        return -1;
    }

    // ------------------------------------
    
    glutInit( &argc, argv );
    glutInitWindowSize( rx, ry );

    glutInitDisplayMode( GLUT_RGBA );
    glutCreateWindow( "physical space viewer" );

    glutReshapeFunc( reshape );
    glutDisplayFunc( display );
    glutMouseFunc( mouse );
    glutMotionFunc( motion );
    glutKeyboardFunc( key );

    init();
    
    if( saveonly )
    {
        reshape( rx, ry );
        display();
        save();
    }
    else
        glutMainLoop();

    return 0;
}