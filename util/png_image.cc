/**
 * No rights reserved, take it if it is useful to you.
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */
#include "png_image.h"

PNGImage::PNGImage(const std::string &img_file, int width, int height, char *title) {
  // Open file for writing (binary mode)
  file_ = fopen(img_file.c_str(), "wb");
  if (file_ == nullptr) {
    fprintf(stderr, "Could not open file %s for writing\n", img_file.c_str());
    return;
  }

  // Initialize write structure
  png_ptr_ = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
  if (png_ptr_ == nullptr) {
    fprintf(stderr, "Could not allocate write struct\n");
    return;
  }

  // Initialize info structure
  info_ptr_ = png_create_info_struct(png_ptr_);
  if (info_ptr_ == nullptr) {
    fprintf(stderr, "Could not allocate info struct\n");
    return;
  }

  // Setup Exception handling
  if (setjmp(png_jmpbuf(png_ptr_))) {
    fprintf(stderr, "Error during png creation\n");
    return;
  }

  png_init_io(png_ptr_, file_);

  // Write header (8 bit colour depth)
  png_set_IHDR(png_ptr_, info_ptr_, width, height,
      8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
      PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

  // Set title
  if (title != nullptr) {
    png_text title_text;
    title_text.compression = PNG_TEXT_COMPRESSION_NONE;
    title_text.key = (char *)"Title";
    title_text.text = title;
    png_set_text(png_ptr_, info_ptr_, &title_text, 1);
  }

  png_write_info(png_ptr_, info_ptr_);

  is_valid_ = true;
}

PNGImage::~PNGImage() {
  Finish();
}

void PNGImage::WriteRow(png_bytep row_buf) {
  png_write_row(png_ptr_, row_buf);
}

void PNGImage::Finish() {
  if (!is_valid_) {
    return;
  }

  if (info_ptr_ != nullptr) {
    png_free_data(png_ptr_, info_ptr_, PNG_FREE_ALL, -1);
    info_ptr_ = nullptr;
  }
  if (png_ptr_ != nullptr) {
    png_write_end(png_ptr_, nullptr);
    png_destroy_write_struct(&png_ptr_, (png_infopp)nullptr);
    png_ptr_ = nullptr;
  }

  if (file_ != nullptr) {
    fclose(file_);
    file_ = nullptr;
  }

  is_valid_ = false;
}
