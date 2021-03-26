// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_FILE_INPUT_H_
#define INCLUDE_FILE_INPUT_H_

/**************************************************************
 * 
 *
 *
 **************************************************************/

#include <string>
#include <fstream>
#include <vector>

std::string find_next_hash(std::ifstream &file_ptr);
std::vector<std::vector<std::string>> read_csv(std::ifstream &file_ptr);
std::vector<std::vector<std::string>> read_ssv(std::ifstream &file_ptr);
std::vector<int> read_itime(std::ifstream &file_ptr, std::string hash);
std::string make_lower(std::string instring);
std::string strip_string_end(std::string instring);
std::string strip_spaces(std::string instring);
std::string read_string(std::ifstream &file_ptr, std::string hash);
int read_int(std::ifstream &file_ptr, std::string hash);
float read_float(std::ifstream &file_ptr, std::string hash);
std::vector<std::string> parse_csv_row_into_vector(std::string line);

#endif  // INCLUDE_FILE_INPUT_H_
