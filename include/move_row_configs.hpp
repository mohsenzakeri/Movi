#ifndef MOVE_ROW_CONFIGS_HPP
#define MOVE_ROW_CONFIGS_HPP

#include <cstdint>

// Large/Constant/Split Index Configuration
#if LARGE_INDEX or CONSTANT_INDEX or SPLIT_INDEX
const uint8_t mask_thresholds1 = static_cast<uint8_t>(~(((1U << 2) - 1) << 0));         // 00000011
const uint8_t mask_thresholds2 = static_cast<uint8_t>(~(((1U << 2) - 1) << 2));         // 00001100
const uint8_t mask_thresholds3 = static_cast<uint8_t>(~(((1U << 2) - 1) << 4));         // 00110000
const uint8_t mask_c = static_cast<uint8_t>(~(((1U << 2) - 1) << 6));                   // 11000000
const uint8_t mask_id =  static_cast<uint8_t>(~(((1U << 4) - 1) << 0));                 // 00001111
const uint8_t mask_overflow_n = static_cast<uint8_t>(~(((1U << 1) - 1) << 4));          // 00010000
const uint8_t mask_overflow_offset = static_cast<uint8_t>(~(((1U << 1) - 1) << 5));     // 00100000
const uint8_t mask_overflow_thresholds = static_cast<uint8_t>(~(((1U << 1) - 1) << 6)); // 01000000
#define MAX_RUN_LENGTH 65535 // 2^16 - 1
#endif

// Regular Index Configuration
#if REGULAR_INDEX
#define SHIFT_C 13
#define SHIFT_ID 12
#define ID_SIG_BITS 4
#define LENGTH_BITS 12
#define C_BITS 3
const uint16_t mask_id =  static_cast<uint16_t>(~(((1U << ID_SIG_BITS) - 1) << SHIFT_ID)); // 11110000 00000000
const uint16_t mask_offset =  static_cast<uint16_t>(~(((1U << LENGTH_BITS) - 1) << 0));    // 00001111 11111111
const uint16_t mask_n =  static_cast<uint16_t>(~(((1U << LENGTH_BITS) - 1) << 0));         // 00001111 11111111
const uint16_t mask_c = static_cast<uint16_t>(~(((1U << C_BITS) - 1) << SHIFT_C));         // 11100000 00000000
#define MAX_RUN_LENGTH 4095 // 2^12-1
#endif

// Regular Thresholds Index Configuration
#if REGULAR_THRESHOLDS_INDEX
#define SHIFT_C 13
#define SHIFT_ID 12
#define SHIFT_THRESHOLD_1 11
#define SHIFT_THRESHOLD_2 11
#define SHIFT_THRESHOLD_3 12
#define ID_SIG_BITS 4
#define LENGTH_BITS 11
#define C_BITS 3
const uint16_t mask_id = static_cast<uint16_t>(~(((1U << ID_SIG_BITS) - 1) << SHIFT_ID));         // 11110000 00000000
const uint16_t mask_offset = static_cast<uint16_t>(~(((1U << LENGTH_BITS) - 1) << 0));            // 00000111 11111111
const uint16_t mask_n = static_cast<uint16_t>(~(((1U << LENGTH_BITS) - 1) << 0));                 // 00000111 11111111
const uint16_t mask_c = static_cast<uint16_t>(~(((1U << C_BITS) - 1) << SHIFT_C));                // 11100000 00000000
const uint16_t mask_thresholds1 = static_cast<uint16_t>(~(((1U << 1) - 1) << SHIFT_THRESHOLD_1)); // 00001000 00000000
const uint16_t mask_thresholds2 = static_cast<uint16_t>(~(((1U << 1) - 1) << SHIFT_THRESHOLD_2)); // 00001000 00000000
const uint16_t mask_thresholds3 = static_cast<uint16_t>(~(((1U << 1) - 1) << SHIFT_THRESHOLD_3)); // 00010000 00000000
#define MAX_RUN_LENGTH 2047 // 2^11-1
#endif

// Blocked Index Configuration
#if BLOCKED_INDEX
#define SHIFT_ID1 10
#define SHIFT_ID2 14
#define SHIFT_ID1_RES 22
#define SHIFT_N 0
#define SHIFT_OFFSET 0
#define SHIFT_C 10
#define ID_SIG_BITS1 6
#define ID_SIG_BITS2 2
#define LENGTH_BITS 10
#define C_BITS 3
const uint16_t mask_id1 = static_cast<uint16_t>(~(((1U << ID_SIG_BITS1) - 1) << SHIFT_ID1));        // 11111100 00000000
const uint16_t mask_id2 = static_cast<uint16_t>(~(((1U << ID_SIG_BITS2) - 1) << SHIFT_ID2));        // 11000000 00000000
const uint16_t mask_offset =  static_cast<uint16_t>(~(((1U << LENGTH_BITS) - 1) << SHIFT_OFFSET));  // 00000011 11111111
const uint16_t mask_n =  static_cast<uint16_t>(~(((1U << LENGTH_BITS) - 1) << SHIFT_N));            // 00000011 11111111
const uint16_t mask_c =  static_cast<uint16_t>(~(((1U << C_BITS) - 1) << SHIFT_C));                 // 00011100 00000000
#define MAX_RUN_LENGTH  1023     // 2^10 - 1
#define BLOCK_SIZE      4194304  // 2^22 -- the actual value might be different
#define MAX_ALLOWED_BLOCKED_ID  16777215 // 2^24 - 1 -- the actual value might be different
#endif

// Blocked Thresholds Index Configuration
#if BLOCKED_THRESHOLDS_INDEX
// Current implementation of MODE 8, does not allocate any bits in the offset for the id field.
// SHIFT_ID2, ID_SIG_BITS2, mask_id2 are only kept for consistency between MODE 2 and MODE 8,
// otherwise they can be removed for MODE 8.
#define SHIFT_ID1 10
#define SHIFT_ID2 16
#define SHIFT_ID1_RES 22
#define SHIFT_OFFSET 0
#define SHIFT_N 0
#define SHIFT_C 10
#define SHIFT_THRESHOLD_1 13
#define SHIFT_THRESHOLD_2 14
#define SHIFT_THRESHOLD_3 15
#define ID_SIG_BITS1 6
#define ID_SIG_BITS2 0
#define LENGTH_BITS 10
#define C_BITS 3
const uint16_t mask_id1 = static_cast<uint16_t>(~(((1U << ID_SIG_BITS1) - 1) << SHIFT_ID1));      // 11111100 00000000
const uint16_t mask_id2 = static_cast<uint16_t>(~(((1U << ID_SIG_BITS2) - 1) << SHIFT_ID2));      // 00000000 00000000
const uint16_t mask_offset = static_cast<uint16_t>(~(((1U << LENGTH_BITS) - 1) << SHIFT_OFFSET)); // 00000011 11111111
const uint16_t mask_n = static_cast<uint16_t>(~(((1U << LENGTH_BITS) - 1) << SHIFT_N));           // 00000011 11111111
const uint16_t mask_c = static_cast<uint16_t>(~(((1U << C_BITS) - 1) << SHIFT_C));                // 00011100 00000000
const uint16_t mask_thresholds1 = static_cast<uint16_t>(~(((1U << 1) - 1) << SHIFT_THRESHOLD_1)); // 00100000 00000000
const uint16_t mask_thresholds2 = static_cast<uint16_t>(~(((1U << 1) - 1) << SHIFT_THRESHOLD_2)); // 01000000 00000000
const uint16_t mask_thresholds3 = static_cast<uint16_t>(~(((1U << 1) - 1) << SHIFT_THRESHOLD_3)); // 10000000 00000000
#define MAX_RUN_LENGTH          1023    // 2^10 - 1
#define BLOCK_SIZE              1048576 // 2^20 (the actual value might be different)
#define MAX_ALLOWED_BLOCKED_ID  4194303 // 2^22 - 1 (the actual value might be different)
#endif

// Tally Index Configuration
#if TALLY_INDEX
#define LENGTH_MASK_BITS 2
#define SHIFT_OFFSET 0
#define SHIFT_N 2
#define C_BITS 4
#define SHIFT_C 4
const uint8_t mask_offset =  static_cast<uint8_t>(~(((1U << LENGTH_MASK_BITS) - 1) << SHIFT_OFFSET));   // 00000011
const uint8_t mask_n =  static_cast<uint8_t>(~(((1U << LENGTH_MASK_BITS) - 1) << SHIFT_N));             // 00001100
const uint8_t mask_c = static_cast<uint8_t>(~(((1U << C_BITS) - 1) << SHIFT_C));                        // 11110000
#define MAX_RUN_LENGTH 1023 // 2^10-1
#endif

// Tally Thresholds Index Configuration
#if TALLY_THRESHOLDS_INDEX
#define LENGTH_MASK_BITS 1
#define SHIFT_OFFSET 0
#define SHIFT_N 1
#define SHIFT_THRESHOLD_1 5
#define SHIFT_THRESHOLD_2 6
#define SHIFT_THRESHOLD_3 7
#define C_BITS 3
#define SHIFT_C 2
const uint8_t mask_offset = static_cast<uint8_t>(~(((1U << LENGTH_MASK_BITS) - 1) << SHIFT_OFFSET)); // 00000001
const uint8_t mask_n = static_cast<uint8_t>(~(((1U << LENGTH_MASK_BITS) - 1) << SHIFT_N));           // 00000010
const uint8_t mask_c = static_cast<uint8_t>(~(((1U << C_BITS) - 1) << SHIFT_C));                     // 00011100
const uint8_t mask_thresholds1 = static_cast<uint8_t>(~(((1U << 1) - 1) << SHIFT_THRESHOLD_1));      // 00100000
const uint8_t mask_thresholds2 = static_cast<uint8_t>(~(((1U << 1) - 1) << SHIFT_THRESHOLD_2));      // 01000000
const uint8_t mask_thresholds3 = static_cast<uint8_t>(~(((1U << 1) - 1) << SHIFT_THRESHOLD_3));      // 10000000
#define MAX_RUN_LENGTH 511 // 2^9-1
#endif

#endif // MOVE_ROW_CONFIGS_HPP