#ifndef FDataSetIO_hh
#define FDataSetIO_hh

// common
#include "FPluginFactory.hh"
#include "FPlugin.hh"
#include "FConfig.hh"

// stl
#include <list>
#include <string>
#include <cassert>
#include <iostream>

// boost
#include <boost/shared_ptr.hpp>

using namespace boost;

// external
class FTensorField;

// forward references
class FTensorFieldIO;
class FTensorFieldInput;
class FTensorFieldOutput;

/**
 * Here we got: this is your interface for writing new
 * input or output plugins. Just inherit from
 * FTensorFieldInput (if you read data into FAnToM)
 * FTensorFieldOutput (if you write data to files)
 *
 * Take a look at dataSet/IO/ for some examples.
 */
class FTensorFieldIO : public FObject
{
    public:
        virtual const std::string getName() const = 0;
};

class FTensorFieldInput : public FTensorFieldIO
{
    public:
        virtual std::list<shared_ptr< FTensorField > > load( const std::list< std::string > &files) =0;
};

class FTensorFieldOutput : public FTensorFieldIO
{
    public:
        virtual bool save( const std::string & file, std::list< shared_ptr< FTensorField > > fields) =0;
};


//! This is just implementation part and should be ignored
//! by the reader *g*


class FIOPlugin;
class FIOPluginFactory;
class FDataSetIO;


/**
 * FIOPlugin is used to call the init function
 * in an object oriented way.
 */
class FIOPlugin : public FObject // later perhaps: FPlugin
{
    public:
        typedef void(*fp_init_fkt)(void);

        FIOPlugin( FDynamicPlugin *p )
        {
            //assert(p->function_map.size() == 1);
            my_init = (fp_init_fkt)(p->getFuncAddress("init"));
        }
        
        FIOPlugin( std::map<std::string, void*> functions)
        {
            //assert(functions.size() == 1);
            my_init = (fp_init_fkt)(functions["init"]);
        }

        void callInit()
        {
            (*my_init)();
        }
    private:
        fp_init_fkt my_init;
};


/**
 * For IO Plugins, the void init( void ) function
 * is called. There every action must be taken
 * by the plugin itself.
 */
class FIOPluginFactory : public FPluginFactory
{
    public:
        FIOPluginFactory() : FPluginFactory()
        { 
        }

        void initialize()
        {
            std::list<std::string> symbols;
            symbols.clear();
            symbols.push_back( "init" );
            
            std::list<std::string> dirs;
            dirs.push_back(FConfig::theConfig()->getEntry("FANTOMBASEDIR")+ "/algoStore/IO/");
            loadPlugins( dirs, symbols );
            // now plugins contains our dynamic plugins
            initPlugins();
        }
        
        // calls init functions of all plugins
        void initPlugins() const
        {
            std::cout << "initializing " << plugins.size() << " IO plugins." << std::endl;
            for(std::vector<FDynamicPlugin*>::const_iterator plugin = plugins.begin(); 
                    plugin != plugins.end(); ++plugin)
            {
                FIOPlugin(*plugin).callInit();
            }
        }
};


/**
 * The new interface for Tensor field input and output operations.
 */
class FDataSetIO : public FObject
{
        typedef shared_ptr<FTensorFieldInput>  rdt;
        typedef shared_ptr<FTensorFieldOutput> wdt;
        typedef std::list<rdt> rlt;
        typedef std::list<wdt> wlt;
        typedef std::list<shared_ptr<FTensorFieldInput> >::iterator rdi;
        typedef std::list<shared_ptr<FTensorFieldOutput> >::iterator wdi;
        
    public:
        static FDataSetIO *theDataSetIO();
    private:
        FDataSetIO();

        void initialize(void)
        {
            pluginFactory.initialize();
        };
    public:        
        ~FDataSetIO();
        
        std::list<std::string> getReaders() const;

        std::list<std::string> getWriter() const;
       
        
        bool isInReaderMap(const std::string & name) const;
        
        bool isInWriterMap(const std::string & name) const;
        
    public: // Plugin Interface, these methods have to be called with appropriate classes
        
        bool registerReader( shared_ptr<FTensorFieldInput> in );

        bool registerWriter( shared_ptr<FTensorFieldOutput> out );
       
    public: // Algorithm interface
        std::list<shared_ptr< FTensorField > > load( const std::string & readername, std::list< std::string > files);
        std::list<shared_ptr< FTensorField > > load( const std::string & readername, std::string files);
        bool save( const std::string & writername, const std::string &files, std::list< shared_ptr<FTensorField> >);
        
    private:

        shared_ptr<FTensorFieldInput> getTensorFieldInput( std::string name ) const;
        shared_ptr<FTensorFieldOutput> getTensorFieldOutput( std::string name ) const;

        FIOPluginFactory pluginFactory;
        rlt  reader;
        wlt writer;
//        mutable bool initialized;
        mutable std::list<std::string > knownPaths;

        static FDataSetIO *dataSetIO;
};

//extern FDataSetIO *theDataSetIO;


#endif
