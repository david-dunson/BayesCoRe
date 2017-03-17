#include <arrayfire.h>
#include <iostream>

using namespace af;
using namespace std;

typedef af::array afarray;

#include "io.hpp"
#include "qnorm.hpp"

#include "rdist_gpu.hpp"
#include "visualize.hpp"

#include "stiefelTensor.hpp"

int main(int argc, char *argv[]) {
  try {
    // Select a device and display arrayfire info
    int device = argc > 1 ? atoi(argv[1]) : 0;
    af::setDevice(device);
    af::info();

    // afarray mu = constant(0, dim4(100, 100));
    // afarray sigma = constant(1, dim4(100, 100));

    // afarray p = randu(dim4(100, 100));
    // afarray z = qnorm(p, mu, sigma);
    // writeCSV(conv_to_std_vector<float>(p), "test/p.csv");
    // writeCSV(conv_to_std_vector<float>(z), "test/z.csv");

    afarray A = af::readArray("../../data/A.af", "A");

    A = A.slices(0, 20);
    A = A.rows(0, 99);
    A = A.cols(0, 99);
    cout << "here" << endl;

    stiefelTensor stTsr(A, 50);
    cout << "here" << endl;

    stTsr.Run(1000);
    writeCSV(conv_to_std_vector<float>(stTsr.D), "test/D.csv");
    writeCSV(conv_to_std_vector<float>(stTsr.Z), "test/Z.csv");

    // const int n1 = 1105;
    // const int n2 = 1105;
    // const int p = 44;
    // const int halfN = n1 / 2;
    // writeCSV(conv_to_std_vector<float>(A.as(f32)), "test/A.csv");

    //  afarray rowIdx = range(dim4(n1, n1));
    // afarray avgA = (sum(A, 2)).as(f32) / (float)A.dims(2);
    // afarray vertex_idx = range(dim4(n1));

    // plotHist(z);

    // plotHist(  z);    //}
  } catch (af::exception &e) {
    fprintf(stderr, "%s\n", e.what());
    throw;
  }

  cout << "Finished" << endl;
  return 0;
}
