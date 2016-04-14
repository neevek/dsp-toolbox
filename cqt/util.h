/**
 * No rights reserved, take it if it is useful to you.
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */
#ifndef UTIL_H_
#define UTIL_H_
#include <png.h>

void val_to_rgb(float level, png_byte *rgb) {
  rgb[0] = rgb[1] = rgb[2] = 255 - fmin(level * 255, 255);
  //rgb[1] = fmin(level * 255, 255);
  //rgb[2] = fmin(level * 255, 255);
} 

#endif /* end of include guard: UTIL_H_ */
