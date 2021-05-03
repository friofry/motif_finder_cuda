#pragma once

#include <cstdint>

uint8_t convert_to_code(char x);
uint8_t to_num(char x);
char to_symbol(uint8_t num);
char hash_letter_to_symbol(uint8_t hash_letter);

uint8_t num_to_hash(uint8_t x);
uint8_t nucl_hash_to_num(uint8_t hash);

