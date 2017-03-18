#include <boost/tokenizer.hpp>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

std::string load_file(const std::string& filename) {
  std::ifstream infile(filename.c_str(), ios::binary);
  std::istreambuf_iterator<char> begin(infile), end;
  return std::string(begin, end);
}

template <class T>
std::vector<T> readCSV(const std::string& filename) {
  std::vector<T> v;

  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

  std::string str = load_file(filename);

  boost::char_separator<char> sep(",\n");
  tokenizer tokens(str, sep);

  for (const auto& t : tokens) {
    // for (tokenizer::iterator tok_iter = tokens.begin();
    // tok_iter != tokens.end(); ++tok_iter){
    v.push_back((T)atof(t.c_str()));
  }

  return v;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& obj) {
  for (auto it = obj.begin(); it != obj.end(); it++) {
    os << *it << endl;
  }
  return os;
}

template <class T>
std::ostream& operator<<(std::ostream& os,
                         const std::vector<std::vector<T>>& obj) {
  for (auto it = obj.begin(); it != obj.end(); it++) {
    std::vector<T> local_obj = *it;
    for (auto it2 = local_obj.begin(); it2 != local_obj.end(); it2++) {
      os << *it2 << ",";
    }
    os << endl;
  }
  return os;
}

template <class T>
std::vector<T> conv_to_std_vector(af::array af_vec) {
  T* af_host_ptr = af_vec.host<T>();
  af::dim4 dim = af_vec.dims();
  int n = dim[0] * dim[1] * dim[2] * dim[3];
  std::vector<T> temp_vector(n);
  temp_vector.assign(af_host_ptr, af_host_ptr + n);
  return temp_vector;
}

template <class T>
void writeCSV(T outputVec, const std::string& filename) {
  std::ofstream of(filename);
  of << outputVec << endl;
  of.close();
}

std::vector<float> operator/(const std::vector<float> obj, const float div) {
  std::vector<float> res(obj.size());
  for (auto i = 0; i < obj.size(); ++i) {
    res[i] = obj[i] / div;
  }
  return res;
}

std::vector<float> operator/(const std::vector<float> obj,
                             const std::vector<float> div) {
  assert(obj.size() == div.size());
  std::vector<float> res(obj.size());
  for (auto i = 0; i < obj.size(); ++i) {
    res[i] = obj[i] / div[i];
  }
  return res;
}
