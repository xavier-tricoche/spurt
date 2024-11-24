#include <iostream>
#include <teem/nrrd.h>
#include <vector>
#include <fstream>

Nrrd* readNrrd(const std::string& filename)
{
    Nrrd *nin = nrrdNew();

    if (nrrdLoad(nin, filename.c_str(), NULL)) {
        char *err = biffGetDone(NRRD);
        std::cerr << "FDTICoherence: " << err << std::endl;
        exit(-1);
    }

    return nin;
}

void bresenham(std::vector< int >& line, int x0, int x1, int y0, int y1)
{
    line.clear();
    int deltax = x1 - x0;
    if (deltax == 0) {
        for (int y = y0 ; y <= y1 ; ++y) {
            line.push_back(x0);
            line.push_back(y);
        }
        return;
    }

    int deltay = y1 - y0;
    float err = 0;
    float deltaerr = deltay / deltax;    // Assume deltax != 0 (line is not vertical),
    int y  = y0;
    for (int x = x0 ; x <= x1 ; ++x) {
        line.push_back(x);
        line.push_back(y);
        err += deltaerr;
        if (fabs(err) >= 0.5) {
            ++y;
            err -= 1.0f;
        }
    }
}

void moments(float& mean, float& var, const std::vector< float >& line, int x, int width)
{
    int min = std::max(x - width, 0);
    int max = std::min(x + width, (int)line.size() - 1);
    
    if (min == max) {
        mean = line[min];
        var = 0;
        mean = 0;
        return;
    }

    mean = 0;
    for (int i = min ; i <= max ; ++i) {
        mean += line[i];
    }
    mean /= (float)(max - min + 1);

    var = 0;
    for (int i = min ; i <= max ; ++i) {
        var += (line[i] - mean) * (line[i] - mean);
    }
    var /= (float)(max - min + 1);
}

int main(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "USAGE: " << argv[0] << " <in> <out> <region width> [<scale>]\n";
        return -1;
    }
    int width = atoi(argv[3]);

    float scale = 0.5;
    if (argc == 5) scale = atof(argv[4]);

    // load image
    Nrrd *nin = readNrrd(argv[1]);
    float *image = (float*)nin->data;
    unsigned int m = nin->axis[0].size;
    unsigned int n = nin->axis[1].size;

    // copy gray scale image into RGB buffer
    float *output = (float*)calloc(m * n * 3, sizeof(float));
    for (unsigned int i = 0 ; i < m*n ; ++i) {
        output[3*i] = output[3*i+1] = output[3*i+2] = image[i];
    }

    // compute max intensity
    float maxcol = *std::max_element(&image[0], &image[m*n]);
    float mincol = *std::min_element(&image[0], &image[m*n]);

    // sample image intensity values along diagonal
    unsigned int N = std::min(m, n);
    float avg = 0;
    std::vector< float > line(N, 0);
    for (unsigned int i = 0 ; i < N ; ++i) {
        unsigned int id = i + i * m; // i-th diagonal point
        line[i] = image[id];
        avg += line[i];
    }
    avg /= (float)N;

    // compute line range for scaling
    float min = *std::min_element(line.begin(), line.end());
    float max = *std::max_element(line.begin(), line.end());

    // draw plot
    int lastx = 0, lasty = 0;
    int x, y;
    std::vector< int > link;
    for (unsigned int i = 0 ; i < N ; ++i) {
        float __min, __max;
        int sign;
        if (line[i] > avg) {
            sign = 1;
        }
        else {
            sign = -1;
        }

        float u = (line[i] - min) / (max - min);
        int shift = floor(u * scale * (float)(m / 2));
        int ii = i + shift;
        // bresenham(link, ii, i, lastx, lasty);
        // for (unsigned int k = 0 ; k < link.size() / 2 ; ++k) {
        // int x = link[2*k];
        // int y = link[2*k+1];
        int x = ii;
        int y = i;
        if (x >= 0 && x < n) {
            int id = x + y * m;
            output[3*id] = maxcol;
            output[3*id+1] = output[3*id+2] = mincol;
        }
        // }
        int id = i + i * m;
        output[3*id] = output[3*id+1] = mincol;
        output[3*id+2] = maxcol;
        lastx = ii;
        lasty = i;
    }

    Nrrd *nout = nrrdNew();
    nrrdWrap_va(nout, output, nrrdTypeFloat, 3, 3, m, n);
    if (nrrdSave(argv[2], nout, NULL)) {
        std::cout << "ERROR while exporting file: " << biffGetDone(NRRD)
        << std::endl;
    }
    nrrdNuke(nout);
    nrrdNuke(nin);

    return 0;
}




































