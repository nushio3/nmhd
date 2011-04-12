#include <fstream>
#include <iostream>
using namespace std;

#include "main_detail.h"
#include "library/cmdline.h"

int main (int argc, char** argv) {
  cmdline::parser p;
  p.add("create-canonical", 'c');
  p.add<string>("output-filename", 'o', "filename to output the benchmark result", false, "");
  p.add<int>("gpu-id", 'g', "gpu id to use", false, 0);

  if(argc <= 1 || !p.parse(argc, argv)) {
    cerr << p.error_full() << p.usage();
    return 1;
  }

  if (p.exist("create-canonical")) {
    main_canonical();
  }
  if (p.exist("output-filename")) {
    main_bench(p.get<int>("gpu-id"), p.get<string>("output-filename"));
  }
  return 0;
}
