#include "toml.h"
#include "cpptoml.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

struct TOMLWrapper {
    std::shared_ptr<cpptoml::table> table;
    TOMLWrapper(std::shared_ptr<cpptoml::table> table) : table(table) {}
};

TOMLWrapper * toml_parse_file(const char * filename) {
    auto config = cpptoml::parse_file(filename);
    return new TOMLWrapper(config);
}

bool toml_contains(TOMLWrapper & wrapper, const char * key) {
    cpptoml::table & table = *(wrapper.table);
    auto exists =  table.contains_qualified(key);
    return exists;
}

#define TOML_GET_TYPE(T) bool toml_get_##T(                  \
    TOMLWrapper & wrapper, const char * key, T & value       \
) {                                                          \
    cpptoml::table & table = *(wrapper.table);               \
    auto option = table.get_qualified_as<T>(key);            \
    if(option) {                                             \
        value = *option;                                     \
    }                                                        \
    return bool(option);                                     \
}
TOML_GET_TYPE(int)
TOML_GET_TYPE(double)
TOML_GET_TYPE(bool)

bool toml_get_string(TOMLWrapper * wrapper, const char * key, char** value, int * length) {
    cpptoml::table & table = *(wrapper->table);
    auto option = table.get_qualified_as<std::string>(key);
    if(option) {
        *length = (*option).length();
        *value = new char[*length + 1];
        strcpy(*value, (*option).c_str());
    }
    return bool(option);
}

void toml_string_delete(char * string) {
    delete[] string;
}
