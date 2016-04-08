#include "CQSpectrogram.h"
#include <stdio.h>
#include <png.h>
#include <iostream>
#include <string>
#include <functional>
#include <algorithm>


class PNGImage {
 public:
   PNGImage(const std::string &img_file, int width, int height, char *title) {
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

   ~PNGImage() {
     if (file_ != nullptr) fclose(file_);
     if (info_ptr_ != nullptr) png_free_data(png_ptr_, info_ptr_, PNG_FREE_ALL, -1);
     if (png_ptr_ != nullptr) png_destroy_write_struct(&png_ptr_, (png_infopp)nullptr);
   }

   bool IsValid() const {
     return is_valid_;
   }

   void WriteRow(png_bytep row_buf) {
     png_write_row(png_ptr_, row_buf);
   }

   void FinishWrite() {
     png_write_info(png_ptr_, info_ptr_);
   }


 private:
   FILE *file_{nullptr};
   png_structp png_ptr_{nullptr};
   png_infop info_ptr_{nullptr};

   bool is_valid_{false};
};

static void val2rgb2(float level, png_byte *rgb) {
  rgb[0] = level * 255;
  rgb[1] = level * 255;
  rgb[2] = level * 255;
}

static void generate_cq_spectrogram() {
  int stereo = 2;

  const int bins_per_octave = 12;

  CQParameters params(44100, 61.735f, 4186.0f, bins_per_octave);
  params.atomHopFactor = 1.f;
  //params.q = 0.3f;
  CQSpectrogram cq(params, CQSpectrogram::InterpolateZeros);

  std::cout << "octaves: " << cq.getOctaves() 
    << ", minFreq=" << cq.getMinFrequency()
    << ", maxFreq=" << cq.getMaxFrequency()
    << ", latency=" << cq.getLatency()
    << ", height=" << (cq.getOctaves() * bins_per_octave)
    << std::endl;

  int frame_size = 1024 * 50;

  short *inbuf = (short *)malloc(sizeof(short) * frame_size * 2); // for stereo
  FILE *file = fopen("/Users/neevek/Desktop/audioprocessing/note_scale/acou/28-acoustic-guitar.raw", "rb");

  CQBase::RealSequence real_seq;
  CQBase::RealBlock real_blk;

  while (1) {
    if (stereo) {
      int n = fread(inbuf, sizeof(short) * 2, frame_size, file);
      if (n != frame_size) {
        break;
      }

      for (int i = 0; i < frame_size; ++i) {
        real_seq.push_back(inbuf[i * 2] + inbuf[i * 2 + 1]);
      //printf("%d\n", inbuf[i]);
      }
    } else {
      int n = fread(inbuf, sizeof(short), frame_size, file);
      if (n != frame_size) {
        break;
      }

      for (int i = 0; i < frame_size; ++i) {
        real_seq.push_back((double)inbuf[i]);
      }
    }


      //printf("===========realseq size: %lu\n", real_seq.size());
    CQBase::RealBlock blk = cq.process(real_seq);
    printf("%lu\n", blk.size());
    real_blk.insert(real_blk.end(), blk.begin(), blk.end());

      //printf("===========blk size: %lu\n", blk.size());
    //for(CQBase::RealColumn &col : blk) {
      //for (double d : col) {
        //printf("%f\n", d);
      //}
      //printf("col size: %lu\n", col.size());
    //}


    //for (const CQBase::ComplexColumn &col : blk) {
    //for (const CQBase::Complex &d : col) {
    //std::cout << d.imag() << ", ";
      //}
      //std::cout << std::endl;
    //}

    real_seq.clear();
  }

  printf("Processing audio data done!\n");

  int width = cq.getOctaves() * bins_per_octave;
  int height = real_blk.size();

  // initialize row buffer
  png_bytep row = (png_bytep) malloc(3 * width * sizeof(png_byte));
  PNGImage png_image("cq_spectrogram.png", width, height, nullptr);

  double nmax = 0;
  double nmin = 0;
  double nrange = 0;
  for(CQBase::RealColumn &col : real_blk) {
    for (double d : col) {
      nmax = std::max(d, nmax);
      nmin = std::min(d, nmin);
    }
  }

  nrange = nmax - nmin;

  for(CQBase::RealColumn &col : real_blk) {
    for (int i = col.size() - 1, j = 0; i >= 0; --i, ++j) {
      val2rgb2((col[i] - nmin)/nrange, row+(j*3) );
    }

    png_image.WriteRow(row);
  }

  png_image.FinishWrite();

  printf("Writing image done!\n");

  fclose(file);
  free(inbuf);
}


int main(int argc, const char *argv[]) {
  generate_cq_spectrogram();

  return 0;
}
