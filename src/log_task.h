#pragma once

#include <Rcpp.h>

// for unix-alike machines only
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
#include <unistd.h>
#include <Rinterface.h>
#endif

void log_task(std::string message);

void log_block_index(int block_index, int n_blocks);