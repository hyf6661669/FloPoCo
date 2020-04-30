#ifndef LOGGING_HPP
#define LOGGING_HPP

#include <iostream>
#include <cstdint>

namespace flopoco{
constexpr const char* remove_trailing_filter(const char* input, const char toRemove)
{
	return (*input == toRemove) ? input+1 : input;
}

constexpr const char* remove_trailing_select(const char* input_int, const char* rec_res, const char toRemove)
{
	return (rec_res == (input_int+1)) ? remove_trailing_filter(input_int, toRemove) :
										rec_res;
}

constexpr const char* remove_trailing_rec(const char* input, const char toRemove);
constexpr const char* remove_trailing_not_null(const char* input, const char toRemove)
{

	return remove_trailing_select(input, remove_trailing_rec(input+1, toRemove), toRemove);
}

constexpr const char* remove_trailing_rec(const char* input, const char toRemove)
{
	return (*input == '\0') ? input : remove_trailing_not_null(input, toRemove);
}

template<size_t N>
constexpr const char* remove_trailing(const char (&input)[N], const char toRemove)
{
	return remove_trailing_rec(static_cast<const char*>(input), toRemove);
}

enum class LogLevel : uint8_t
{
	List = 0,
	Info = 1,
	Detailed = 2,
	Debug = 3,
	Full = 4
};

extern int MAXLOGLEVEL;
}

// Reporting levels
#define LIST 0       // information necessary to the user of FloPoCo
#define INFO 1       // information useful to the user of FloPoCo
#define DETAILED 2   // information that shows how the algorithm works, useful to inquisitive users
#define DEBUG 3      // debug info, useful mostly to developers
#define FULL 4       // pure noise
#define REPORT(level, stream) {if ((level)<=MAXLOGLEVEL){ cerr << "> " << std::string(remove_trailing(__FILE__, '/')) << " " << uniqueName_ <<": " << stream << endl;}else{}} 
#define THROWERROR(stream) {{ostringstream o; o << " ERROR in " << uniqueName_ << " (" << srcFileName << "): " << stream << endl; throw o.str();}}
#endif
