#pragma once

#include <string>
#include <cstdint>


enum ErrorCode { OK, LICENSE_ERROR };
static const uint16_t LICENSE_DATA[3] = { 2026^0x55, 2^0x55, 1^0x55 };
ErrorCode check_license();
std::string get_license_info();

static void decode_license(int &year, int &month, int &day) {
    year  = LICENSE_DATA[0] ^ 0x55;
    month = LICENSE_DATA[1] ^ 0x55;
    day   = LICENSE_DATA[2] ^ 0x55;
}
