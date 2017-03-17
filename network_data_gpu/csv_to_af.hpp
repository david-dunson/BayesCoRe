//#include "io.hpp"

void csv_to_af() {
  const int p = 44;

  std::vector<int> A_host = readCSV<int>("BNU1.csv");

  std::cout << "Finish Reading" << std::endl;
  af::array A(dim4(n1, n2, p), &A_host[0], afHost);

  af::saveArray("A", A, "A.af");
}