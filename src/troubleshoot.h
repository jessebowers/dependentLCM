#pragma once

// Are we troubleshooting? 0 for no. 1 for print. 2 for runtimes.
#define TROUBLESHOOT 0

// only compile if troubleshooting...
#if TROUBLESHOOT > 0

#include <chrono>
typedef std::chrono::time_point<std::chrono::high_resolution_clock> ttime;

// functions
void trouble_start_print(const std::string& function_name);
void trouble_end_print(const std::string& function_name);
void trouble_init_print();
Rcpp::List trouble_list_print();
ttime trouble_start_time(const std::string& function_name);
void trouble_end_time(ttime _trouble_start_time, const std::string& function_name);
void trouble_init_time();
Rcpp::List trouble_list_time();

#endif

//***
//*** Define Macros:
//***

#if TROUBLESHOOT == 1

#define TROUBLE_START(str) std::string trouble_function_name = str; trouble_start_print(trouble_function_name);
#define TROUBLE_END ; trouble_end_print(trouble_function_name);
#define TROUBLE_INIT trouble_init_print();
#define TROUBLE_LIST trouble_list_print()

#elif TROUBLESHOOT == 2

#define TROUBLE_START(str) std::string trouble_function_name = str; ttime _trouble_start_time = trouble_start_time(trouble_function_name);
#define TROUBLE_END ; trouble_end_time(_trouble_start_time, trouble_function_name);
#define TROUBLE_INIT trouble_init_time();
#define TROUBLE_LIST trouble_list_time()

#else
// else TROUBLSHOOT=0 and ignore these macros
#define TROUBLE_START(str)
#define TROUBLE_END
#define TROUBLE_INIT
#define TROUBLE_LIST Rcpp::List::create()

#endif