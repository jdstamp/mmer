#include "log_task.h"

void log_task(std::string message) {

        std::cout << "\r" << message << "\r" << std::flush;

        // for unix-alike machines only
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
        R_FlushConsole();
#endif
        R_CheckUserInterrupt();

}

void log_block_index(int block_index, int n_blocks) {
        std::string message = "block " + std::to_string(block_index + 1) + " of " + std::to_string(n_blocks);
        log_task(message);
}