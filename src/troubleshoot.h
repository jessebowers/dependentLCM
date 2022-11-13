#pragma once

// Are we troubleshooting? 0 for no. 1 for print. 2 for runtimes.
#define TROUBLESHOOT 0

// only compile if troubleshooting...
#if TROUBLESHOOT > 0

// functions
unsigned long long int trouble_start_print(std::string function_name);
void trouble_end_print(unsigned long long int trouble_id, std::string function_name);
void trouble_init_print();
Rcpp::List trouble_list_print();
unsigned long long int trouble_start_time(std::string function_name);
void trouble_end_time(unsigned long long int trouble_id, std::string function_name);
void trouble_init_time();
Rcpp::List trouble_list_time();
#endif

//***
//*** Define Macros:
//***

#if TROUBLESHOOT == 1

#define TROUBLE_START(str) std::string trouble_function_name = str; unsigned long long int trouble_id = trouble_start_print(trouble_function_name);
#define TROUBLE_END ; trouble_end_print(trouble_id, trouble_function_name);
#define TROUBLE_INIT trouble_init_print();
#define TROUBLE_LIST trouble_list_print()

#elif TROUBLESHOOT == 2

#define TROUBLE_START(str) std::string trouble_function_name = str; unsigned long long int trouble_id = trouble_start_time(trouble_function_name);
#define TROUBLE_END ; trouble_end_time(trouble_id, trouble_function_name);
#define TROUBLE_INIT trouble_init_time();
#define TROUBLE_LIST trouble_list_time()

#else
// else TROUBLSHOOT=0 and ignore these macros
#define TROUBLE_START(str)
#define TROUBLE_END
#define TROUBLE_INIT
#define TROUBLE_LIST Rcpp::List::create()

#endif