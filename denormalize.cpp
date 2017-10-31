#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <cstring>

int main(int argc, char *argv[]){//data_count x 2
  if(argc != 3){
    std::cout << "usage: " << argv[0] << " normalized.bedGraph counts.bedGraph";
    return 1;
  }
  std::ifstream bedGraph_file(argv[1]);
  if(!bedGraph_file.is_open()){
    std::cout << "Could not open data file\n";
    return 2;
  }
  std::string line;
  int chromStart, chromEnd, items, line_i=0, n_zeros=0;
  char chrom[100];
  char extra[100] = "";
  double coverage, min_coverage = INFINITY;
  while(std::getline(bedGraph_file, line)){
    line_i++;
    items = sscanf
      (line.c_str(),
       "%s %d %d %lf%s\n",
       chrom, &chromStart, &chromEnd, &coverage, extra);
    //printf("%s %d %d %d%s\n", chrom, chromStart, chromEnd, coverage, extra);
    if(items < 4){
      printf("error: expected '%%s %%d %%d %%f\\n' on line %d\n", line_i);
      std::cout << line << "\n";
      return 3;
    }
    if(0 < strlen(extra)){
      printf("error: non-numeric data on line %d\n", line_i);
      std::cout << line << "\n";
      return 4;
    }
    if(0 == coverage){
      n_zeros++;
    }else{
      if(coverage < 0){
	printf("error: negative coverage on line %d\n", line_i);
	std::cout << line << "\n";
	return 5;
      }else{
	if(coverage < min_coverage){
	  min_coverage = coverage;
	}
      }
    }
  }
  bedGraph_file.clear();
  bedGraph_file.seekg(0, std::ios::beg);
  printf("min_coverage=%f n_zeros=%d/%d\n", min_coverage, n_zeros, line_i);
  std::ofstream out_file;
  out_file.open(argv[2]);
  if(!out_file.is_open()){
    std::cout << "Could not open out file\n";
    return 6;
  }
  while(std::getline(bedGraph_file, line)){
    line_i++;
    items = sscanf
      (line.c_str(),
       "%s %d %d %lf\n",
       chrom, &chromStart, &chromEnd, &coverage);
    out_file << chrom <<
      "\t" << chromStart <<
      "\t" << chromEnd <<
      "\t" << coverage/min_coverage <<
      "\n";
  }
  return 0;
}

