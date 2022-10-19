#ifndef FPluginFactory_hh
#define FPluginFactory_hh

// common
#include "FObject.hh"
#include "FPlugin.hh"

// stl
#include <list>
#include <string>
#include <vector>

class FPluginFactory : public FObject
{
    public:
        FPluginFactory();
        bool loadPlugins(std::list<std::string> dirs, std::list<std::string> symbols);
        ~FPluginFactory();

        const FString& getClassName() const;
    protected:
        std::vector<FDynamicPlugin*> plugins;
};

#endif

