#include "exists.hpp"

namespace gimbricate {

bool file_exists(const std::string& name) {
    struct stat buffer;   
    return (stat (name.c_str(), &buffer) == 0); 
}

}
