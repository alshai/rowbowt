#ifndef UTILS_RI_HPP_
#define UTILS_RI_HPP_

#include <sstream>
#include <iostream>
#include <chrono>
#include <cstring>
#include <sys/stat.h>

using namespace std;


namespace ri {
namespace utils {

inline void die(std::string msg) {
    cerr << msg << endl;
    exit(1);
}

inline uint8_t bitsize(uint64_t x){
	if(x==0) return 1;
	return 64 - __builtin_clzll(x);
}

class Timer {

    public:

	Timer(std::string &&name=""): name_{std::move(name)}, 
                                  start_(std::chrono::system_clock::now()) 
                                  {}

	void stop() {
        stop_ = std::chrono::system_clock::now();
    }

	void restart() {
        start_ = std::chrono::system_clock::now();
    }

	void report() {
        std::cerr << name_ << " " << std::chrono::duration_cast<std::chrono::microseconds>(stop_ - start_).count() << "\n";
    }

	~Timer() {
        stop(); 
        report();
    }

	void rename(const char *name) {
        name_ = name;
    }

    private: 
	using TpType = typename std::chrono::system_clock::time_point;
	std::string name_;
	TpType start_, stop_;
};

inline bool file_exists(const std::string& name) {
  struct stat buffer;   
  return (stat(name.data(), &buffer) == 0); 
}

static size_t get_file_size(const char* path) {
    if (!strcmp(path, "-")) {
        fprintf(stderr, "cannot get file size from stdin");
        exit(1);
    }
    FILE* fp = fopen(path, "rb");
    if (fp == NULL) {
        fprintf(stderr, "error opening file");
        exit(1);
    }
    fseek(fp, 0, SEEK_END);
    size_t size = ftell(fp);
    fclose(fp);
    return size;
}

}
}


#endif /* UTILS_RI_HPP_ */
