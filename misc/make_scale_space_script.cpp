#include <iostream>
#include <fstream>
#include <teem/hest_helper.hpp>

char* _base, *output;
int nscales, ninterp;
float minscale, maxscale;

void initialize(int argc, char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "b",      "_base",            airTypeString,  1, 1, &_base,   NULL,       "base name of NRRD scalar field");
    hestOptAdd(&hopt, "o",      "output",               airTypeString,  1, 1, &output,      NULL,       "script output name");
    hestOptAdd(&hopt, "nc",     "# computed scales",    airTypeInt,     0, 1, &nscales,     "15",       "number of computed scales");
    hestOptAdd(&hopt, "ni",     "# interp. scales",     airTypeInt,     0, 1, &ninterp,     "40",       "number of interpolated scales");
    hestOptAdd(&hopt, "min",    "min scale",            airTypeFloat,   0, 1, &minscale,    "0",        "finest scale");
    hestOptAdd(&hopt, "max",    "max scale",            airTypeFloat,   0, 1, &maxscale,    "5",        "coarsest scale");
    
    __hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                     me, "Generate script for scale space interpolation using teem",
                     AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    std::fstream file(output, std::ios::out);
    
    float dscale = (float)(maxscale - minscale)/(float)ninterp;
    
    file << "#!/bin/bash\n";
    // pre-sample scale space
    file << "vprobe -i " << _base << ".nrrd -k scalar -q v -ssn " << nscales << " -ssr " << minscale << " " << maxscale
         << " -ssw 0 -sssf " << _base << "-sigma=\%u.nrrd -ssnd -o " << _base << "-sigma=0-v.nrrd\n";
    // compute value
    for (int i=1 ; i<=ninterp ; ++i) {
        float s = minscale + (float)i*dscale;
        file << "vprobe -i " << _base << ".nrrd -k scalar -q v -ssn " << nscales << " -ssr " << minscale << " " << maxscale
             << " -ssw " << s << " -ssrf " << _base << "-sigma=\%u.nrrd -ssnd -o " << _base << "-sigma=" << i << "-v.nrrd\n";
    }
    // compute gradient
    for (int i=0 ; i<=ninterp ; ++i) {
        float s = minscale + (float)i*dscale;
        file << "vprobe -i " << _base << ".nrrd -k scalar -q gv -ssn " << nscales << " -ssr " << minscale << " " << maxscale
             << " -ssw " << s << " -ssrf " << _base << "-sigma=\%u.nrrd -ssnd -o " << _base << "-sigma=" << i << "-gv.nrrd\n";
    }
    // compute hessian
    for (int i=0 ; i<=ninterp ; ++i) {
        float s = minscale + (float)i*dscale;
        file << "vprobe -i " << _base << ".nrrd -k scalar -q hess -ssn " << nscales << " -ssr " << minscale << " " << maxscale
             << " -ssw " << s << " -ssrf " << _base << "-sigma=\%u.nrrd -ssnd -o " << _base << "-sigma=" << i << "-hess.nrrd\n";
    }
    // compute ridge strength
    for (int i=0 ; i<=ninterp ; ++i) {
        float s = minscale + (float)i*dscale;
        file << "vprobe -i " << _base << ".nrrd -k scalar -q hesseval2 -ssn " << nscales << " -ssr " << minscale << " " << maxscale
             << " -ssw " << s << " -ssrf " << _base << "-sigma=\%u.nrrd -ssnd -o " << _base << "-sigma=" << i << "-hesseval2.nrrd\n";
    }
    // filter ridge strength
    file << "for ((n=0 ; n<=" << ninterp << "; ++n)); do\n"
         << "unu 2op x -1 " << _base << "-sigma=${n}-hesseval2.nrrd | unu 2op max - 0 -o " << _base << "-sigma=${n}-lambda2-r.nrrd\n";
    file << "done\n";
    
    // compute optimal scale
    
    // initialize max volume to zero
    file << "unu 2op x 0 " << _base << "-sigma=0-lambda2-r.nrrd -o " << _base << "-lambda2-max-r.nrrd\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "unu 2op max " << _base << "-lambda2-max-r.nrrd " << _base << "-sigma=$a-lambda2-r.nrrd -o " << _base << "-lambda2-max-r.nrrd\n";
    file << "done\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "unu 2op eq " << _base << "-sigma=$a-lambda2-r.nrrd " << _base << "-lambda2-max-r.nrrd -o mask-$a-r.nrrd\n";
    file << "done\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "unu 2op x mask-$a-r.nrrd $a -o mask-$a-r.nrrd\ndone\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "cp mask-$a-r.nrrd -o mask2-$a-r.nrrd\ndone\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "unu 2op + mask-$a-r.nrrd mask2-$a-r.nrrd -o mask-$a-r.nrrd\ndone\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "unu 2op eq mask-$a-r.nrrd 0 -o test$a.nrrd\ndone\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "unu 2op x 1000 test$a.nrrd -o test$a.nrrd\ndone\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "unu 2op + mask-$a-r.nrrd test$a.nrrd -o mask-$a-r.nrrd\ndone\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "removecontent mask-$a-r.nrrd\ndone\n";
    file << "unu 2op min 1000 mask-40-r.nrrd -o " << _base << "-lambda2-bestscale-r.nrrd\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "unu 2op min mask-$a-r.nrrd " << _base << "-lambda2-bestscale-r.nrrd -o " << _base << "-lambda2-bestscale-r.nrrd\ndone\n";
    file << "unu 2op - " << _base << "-lambda2-bestscale-r.nrrd 1 -o " << _base << "-lambda2-bestscale-r.nrrd\n";
    file << "unu save -f nrrd -e ascii -i " << _base << "-lambda2-bestscale-r.nrrd -o " << _base << "-lambda2-bestscale-r-ascii.nrrd\n";
    file << "unu join -a 0 -incr -i " << _base << "-lambda2-bestscale-r.nrrd " << _base << "-lambda2-bestscale-r.nrrd " << _base << "-lambda2-bestscale-r.nrrd -o " << _base << "-lambda2-bestscale-r-x3.nrrd\n";
    file << "unu join -a 0 -i " << _base << "-lambda2-bestscale-r-x3.nrrd " << _base << "-lambda2-bestscale-r-x3.nrrd " << _base << "-lambda2-bestscale-r-x3.nrrd -o " << _base << "-lambda2-bestscale-r-x9.nrrd\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "unu 2op eq "<< _base << "-lambda2-bestscale-r.nrrd $a -o at$a.nrrd\ndone\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "unu 2op x " << _base << "-sigma=$a-r.nrrd at$a.nrrd -o at$a.nrrd\ndone\n";
    file << "unu 2op x 0 " << _base << "-lambda2-bestscale-r.nrrd -o " << _base << "-val-lambda2-max-r.nrrd\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "unu 2op + " << _base << "-val-lambda2-max-r.nrrd at$a.nrrd -o " << _base << "-val-lambda2-max-r.nrrd\ndone\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "unu 2op eq " << _base << "-lambda2-bestscale-r-x3.nrrd $a -o at$a.nrrd\ndone\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "unu 2op x " << _base << "-sigma=$a-gv-r.nrrd at$a.nrrd -o at$a.nrrd\ndone\n";
    file << "unu 2op x 0 " << _base << "-lambda2-bestscale-r-x3.nrrd -o " << _base << "-gv-lambda2-max-r.nrrd\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "unu 2op + " << _base << "-gv-lambda2-max-r.nrrd at$a.nrrd -o " << _base << "-gv-lambda2-max-r.nrrd\ndone\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "unu 2op eq " << _base << "-lambda2-bestscale-r-x9.nrrd $a -o at$a.nrrd\ndone\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "unu 2op x " << _base << "-sigma=$a-hess-r.nrrd at$a.nrrd -o at$a.nrrd\ndone\n";
    file << "unu 2op x 0 " << _base << "-lambda2-bestscale-r-x9.nrrd -o " << _base << "-hess-lambda2-max-r.nrrd\ndone\n";
    file << "for ((  a = 0 ;  a <= " << ninterp << ";  a++  )); do\n";
    file << "unu 2op + " << _base << "-hess-lambda2-max-r.nrrd at$a.nrrd -o " << _base << "-hess-lambda2-max-r.nrrd\ndone\n";
    file << "unu save -f nrrd -i " << _base << "-lambda2-bestscale-r.nrrd -o scale.nrrd\n";
    file << "unu save -f nrrd -i " << _base << "-val-lambda2-max-r.nrrd -o " << _base << ".nrrd\n";
    file << "unu save -f nrrd -i " << _base << "-gv-lambda2-max-r.nrrd -o grad.nrrd\n";
    file << "unu save -f nrrd -i " << _base << "-hess-lambda2-max-r.nrrd -o hess.nrrd\n";
    file << "unu save -f nrrd -i " << _base << "-lambda2-max-r.nrrd -o evl2.nrrd\n";
    file << "unu 2op x -1 evl2.nrrd -o evl2.nrrd\n";
    file << "rm mask*\n";
    file << "rm test*\n";
    file << "rm at*\n";
    file << "rm tmp*\n";
    
    file.close();
    
    return 0;
}
