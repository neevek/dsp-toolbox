/**
 * No rights reserved, take it if it is useful to you.
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */
#ifndef PNG_IMAGE_H_
#define PNG_IMAGE_H_
#include <string>
#include <png.h>

class PNGImage {
 public:
   PNGImage(const std::string &img_file, int width, int height, char *title);
   ~PNGImage();

   bool IsValid() const {
     return is_valid_;
   }

   void WriteRow(png_bytep row_buf);
   void Finish();

 private:
   FILE *file_{nullptr};
   png_structp png_ptr_{nullptr};
   png_infop info_ptr_{nullptr};

   bool is_valid_{false};
};


#endif /* end of include guard: PNG_IMAGE_H_ */
