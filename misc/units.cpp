#include <iostream>
#include <map>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/cmath.hpp>

// quantities
#include <boost/units/systems/si/area.hpp>
#include <boost/units/systems/si/io.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/mass.hpp>
#include <boost/units/systems/si/temperature.hpp>
#include <boost/units/systems/si/velocity.hpp>
#include <boost/units/systems/si/volume.hpp>

// metric system
#include <boost/units/base_units/metric/liter.hpp>
#include <boost/units/base_units/metric/hectare.hpp>
#include <boost/units/base_units/metric/foot.hpp>
#include <boost/units/base_units/metric/nautical_mile.hpp>
#include <boost/units/base_units/metric/ton.hpp>

// SI units
#include <boost/units/base_units/si/kelvin.hpp>
#include <boost/units/base_units/si/kilogram.hpp>
#include <boost/units/base_units/si/meter.hpp>

// temperature
#include <boost/units/base_units/temperature/celsius.hpp>
#include <boost/units/base_units/temperature/conversion.hpp>
#include <boost/units/base_units/temperature/fahrenheit.hpp>

// US customary units
#include <boost/units/base_units/us/cup.hpp>
#include <boost/units/base_units/us/fluid_once.hpp>
#include <boost/units/base_units/us/foot.hpp>
#include <boost/units/base_units/us/gallon.hpp>
#include <boost/units/base_units/us/inch.hpp>
#include <boost/units/base_units/us/mile.hpp>
#include <boost/units/base_units/us/ounce.hpp>
#include <boost/units/base_units/us/pound.hpp>
#include <boost/units/base_units/us/quart.hpp>
#include <boost/units/base_units/us/tablespoon.hpp>
#include <boost/units/base_units/us/teaspoon.hpp>
#include <boost/units/base_units/us/ton.hpp>
#include <boost/units/base_units/us/yard.hpp>

using namespace boost::units;
std::map<std::string, quantity<si::length> >      known_length_units;
std::map<std::string, quantity<si::area> >        known_area_units;
std::map<std::string, quantity<si::volume> >      known_volume_units;
std::map<std::string, quantity<si::temperature> > known_temperature_units;

void displayUsageAndExit(const std::string& what) {
    if (what.size()) std::cerr << "ERROR: " << what << '\n';
    std::cerr 
        << "USAGE: " << me << " [options]\n"
        << "DESCRIPTION: Convert value from one unit to another\n"
        << "OPTIONS:\n"
        << " -h | --help                  Print this information\n"
        << " -i | --input <string>        Input unit\n"
        << " -o | --output <string>       Output unit\n";
    exit(1);
}
