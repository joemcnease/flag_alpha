#include <ctime>

#include "license.hpp"


ErrorCode check_license() {
    int y, m, d;
    decode_license(y, m, d);

    std::tm expiry{};
    expiry.tm_year = y - 1900;
    expiry.tm_mon  = m - 1;
    expiry.tm_mday = d;
    std::time_t expiry_time = std::mktime(&expiry);

    std::time_t now = std::time(nullptr);
    return now > expiry_time ? LICENSE_ERROR : OK;
}

std::string get_license_info() {
    int y, m, d;
    decode_license(y, m, d);

    char buf[32];
    std::tm expiry{};
    expiry.tm_year = y - 1900;
    expiry.tm_mon  = m - 1;
    expiry.tm_mday = d;
    std::mktime(&expiry);

    std::strftime(buf, sizeof(buf), "%Y-%m-%d", &expiry);
    return buf;
}


