#pragma once

#include "BooPHF.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <random>
#include <algorithm>
#include <fstream>
#include "threads.hpp"

namespace gimbricate {

//example with user provided custom hasher for uint64_t type :

class Custom_string_Hasher
{
public:
	// the class should have operator () with this signature :
	uint64_t operator ()   (std::string key, uint64_t seed=0) const
	{
		uint64_t hash  =  hash_fn(key);
		hash ^= seed;
		return hash;
	}
    std::hash<std::string> hash_fn;
};

//then tell BBhash to use this custom hash : (also appears below, line 104)
typedef boomphf::mphf<  std::string, Custom_string_Hasher  > boophf_t;

class node_index {
    boophf_t* bphf = nullptr;
public:
    node_index(const std::vector<std::string>& input);
    ~node_index(void);
    uint64_t get_id(const std::string& key);
};

}
