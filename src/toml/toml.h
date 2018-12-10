#ifndef _TOML_H_
#define _TOML_H_

class TOMLWrapper;

extern "C" {
    TOMLWrapper * toml_parse_file(const char * filename);

    bool toml_contains(TOMLWrapper & wrapper, const char * key);

    bool toml_get_double(TOMLWrapper & wrapper, const char * key, double & value);
    bool toml_get_int(TOMLWrapper & wrapper, const char * key, int & value);
    bool toml_get_bool(TOMLWrapper & wrapper, const char * key, bool & value);

    bool toml_get_string(TOMLWrapper * wrapper, const char * key, char** value, int * length);
    void toml_string_delete(char * string);
}

#endif
