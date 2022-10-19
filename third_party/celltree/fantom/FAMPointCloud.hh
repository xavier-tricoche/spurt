#ifndef FAMPointCloud_HH
#define FAMPointCloud_HH

// math
#include "FFixArray.hh"

// analysis
#include "FAMElement.hh"

// stl
#include <iostream>
#include <fstream>

class FAMPointCloud : public FAMElement
{
    public:
        FAMPointCloud(){}

        virtual const FString& getClassName() const 
        {
            static FString classname("FAMPointCloud");
            return classname;
        }

        void clear()
        {
            points.clear();
        }

        bool isEmpty() const
        {
            return points.size() == 0;
        }

        const std::vector<FArray3>& getPoints() const
        {
            return points;
        }

        std::vector<FArray3>& getPoints()
        {
            return points;
        }

        void setPoints(std::vector<FArray3>& npoints)
        {
            points =npoints;
        }

        size_t size() const
        {
            return points.size();
        }

        void save(const FString& fileNameBase)
        {
            std::ofstream out( (fileNameBase + ".pc").c_str() );
            out << "[" << getClassName() << "]" << std::endl;
            out << points.size() << std::endl;
            std::copy( points.begin(), points.end(), ostream_iterator<FArray3>(out) );
            out.flush();
            out.close();
        }

        void load(const FString& fileNameBase)
        {
            std::ifstream in( (fileNameBase + ".pc").c_str());
            std::string tmp;
            in >> tmp;
            if(tmp != (std::string("[") + getClassName() + "]"))
                THROW_EXCEPTION( FException, "Invalid Header, wrong file?" );

            this->clear();
            unsigned int size;
            in >> size;
            std::copy( istream_iterator<FArray3>(in), istream_iterator<FArray3>(), back_inserter(points) );
            if(points.size() != size)
                THROW_EXCEPTION( FException, "Invalid File, file damaged?" );
            in.close();
        }

    private:
        std::vector<FArray3> points;
};

#endif
