#include "letter_conversions.h"

#include <vector>

#include <config.h>

using namespace std;
namespace  {

const char _symbols[ALPH_SIZE+1] = {
        //      0xCGTA
        '?',
        'A', //  1 - 0x0001 T
        'T', //  2 - 0x0010 A
        'W', //  3 - 0x0011 W
        'G', //  4 - 0x0100 C
        'R', //  5 - 0x0101 Y
        'K', //  6 - 0x0110 M
        'D', //  7 - 0x0111 H
        'C', //  8 - 0x1000 G
        'M', //  9 - 0x1001 K
        'Y', // 10 - 0x1010 R
        'H', // 11 - 0x1011 D
        'S', // 12 - 0x1100 S
        'V', // 13 - 0x1101 B
        'B', // 14 - 0x1110 V
        'N'  // 15 - 0x1111 N
};
} // namespace

inline uint8_t convert_to_code(char x)
{
	switch(x) {
		case 'A' : return 1;  // 0x0001
		case 'T' : return 2;  // 0x0010
		case 'G' : return 4;  // 0x0100
		case 'C' : return 8;  // 0x1000
		case 'R' : return 5;  // 0x0101
		case 'Y' : return 10; // 0x1010
		case 'M' : return 9;  // 0x1001
		case 'K' : return 6;  // 0x0110
		case 'W' : return 3;  // 0x0011
		case 'S' : return 12; // 0x1100
		case 'B' : return 14; // 0x1110
		case 'H' : return 11; // 0x1011
		case 'V' : return 13; // 0x1101
		case 'D' : return 7;  // 0x0111
		case 'N' : return 15; // 0x1111

		case 'a' : return 1;  // 0x0001
		case 't' : return 2;  // 0x0010
		case 'g' : return 4;  // 0x0100
		case 'c' : return 8;  // 0x1000
		case 'r' : return 5;  // 0x0101
		case 'y' : return 10; // 0x1010
		case 'm' : return 9;  // 0x1001
		case 'k' : return 6;  // 0x0110
		case 'w' : return 3;  // 0x0011
		case 's' : return 12; // 0x1100
		case 'b' : return 14; // 0x1110
		case 'h' : return 11; // 0x1011
		case 'v' : return 13; // 0x1101
		case 'd' : return 7;  // 0x0111
		case 'n' : return 15; // 0x1111
	};
	return 0;
}

uint8_t to_num(char x)
{
	switch(x) {
		case 'a': return 0;
		case 'A': return 0;
		case 't': return 1;
		case 'T': return 1;
		case 'g': return 2;
		case 'G': return 2;
		case 'c': return 3;
		case 'C': return 3;
	};

	return 255;
}

uint8_t num_to_hash(uint8_t x)
{
    return 1 << x;
}

uint8_t nucl_hash_to_num(uint8_t hash)
{
    static const vector<uint8_t> hash_to_num = {
        255, // 0000 0
        0,   // 0001 1
        1,   // 0010 2
        255, // 0011 3
        2,   // 0100 4
        255, // 0101 5
        255, // 0110 6
        255, // 0111 7
        3,   // 1000 8
        255, // 1001 9
        255, // 1010 10
        255, // 1011 11
        255, // 1100 12
        255, // 1101 13
        255, // 1110 14
        255, // 1111 15
    };

    return hash_to_num[hash];
}

char to_symbol(uint8_t num)
{
    switch (num) {
        case 0: return 'A';
        case 1: return 'T';
        case 2: return 'G';
        case 3: return 'C';
    }
    return '-';
}

char hash_letter_to_symbol(uint8_t hash_letter) {
    return _symbols[hash_letter];
}

uint8_t to_compl_num(uint8_t num)
{
    static std::vector<uint8_t> complmap = {1, 0, 3, 2};
    return complmap.at(num);
}
