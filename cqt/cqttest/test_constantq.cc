#include <stdio.h>
#include <string>
#include <png.h>
#include "ConstantQTransform.hxx"
#include "FourierTransform.hxx"

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

// http://en.wikipedia.org/wiki/Window_function
float hamming_window(int n, int N) {
  return 0.53836f - 0.46f * (float)cos((2 * M_PI * n) / (N - 1));
}

static void do_cqt(const char *raw_pcm_file,
                   const char *cqt_spectrogram_filepath,
                   const char *cqt_data_filepath,
                   const char *feature) {
  bool is_stereo = true;

  unsigned int sample_rate = 16000;
  double min_freq = 61.735f;
  double max_freq = 4186.0f;
  unsigned int bins_per_octave = 12;

  double sparse_cqt_kernel_threshold = 0.0054;
  

  Simac::ConstantQTransform cqt(sample_rate, min_freq, max_freq, bins_per_octave);
  unsigned int frame_size = cqt.getfftlength();
  printf("initializing fft...\n");
  FourierTransform fft(frame_size, 1, 0);
  printf("initializing cqt...\n");
  cqt.sparsekernel(sparse_cqt_kernel_threshold);

  printf("frame_size=%d\n", frame_size);

  FILE *file = fopen(raw_pcm_file, "rb");
  short inbuf[frame_size * 2];
  double frame_data[frame_size];

  std::vector<std::vector<double>> cqt_data;

  printf("doing fft & cqt\n");
  while(1) {
    if (is_stereo) {
      int n = fread(inbuf, sizeof(short) * 2, frame_size, file);
      if (n != frame_size) {
        break;
      }

      for (int i = 0; i < frame_size; ++i) {
        //frame_data[i] = inbuf[i * 2] + inbuf[i * 2 + 1];
        // use right channel only in this case
        printf("i=%d, hamming: %f\n", i, hamming_window(i, frame_size));
        frame_data[i] = inbuf[i * 2 + 1] * hamming_window(i, frame_size);
      }

      long overlap_size = (long)(frame_size * sizeof(short) * 2 * 0.5f);
      if (overlap_size % 2 != 0) {
        --overlap_size;
      }
      fseek(file, -overlap_size, SEEK_CUR);
      //fseek(file, ftell(file)-overlap_size, SEEK_SET);
      printf("overlap=%ld, frame_size=%d, cur_pos=%ld\n", overlap_size, frame_size, ftell(file));
    } else {
      int n = fread(inbuf, sizeof(short), frame_size, file);
      if (n != frame_size) {
        break;
      }

      for (int i = 0; i < frame_size; ++i) {
        frame_data[i] = inbuf[i];
      }
    }



    fft.doIt(frame_data);
    cqt.doIt(fft.spectrum());

    std::vector<double> spectrum = cqt.constantQSpectrum();
    std::vector<double> col(spectrum);

    for (int i = 0; i < col.size()/2; ++i) {
      col[i] = sqrt(col[i*2] * col[i*2] + col[i*2+1] * col[i*2+1]);
    }
    col.resize(col.size()/2);

    cqt_data.push_back(col);

    //printf("================================================\n");
    //for (int i = 0; i < spectrum.size(); ++i) {
      //printf("%lf\n", spectrum[i]);
    //}
  }

  if (cqt_data.size() > 0) {
    FILE *fdata = fopen(cqt_data_filepath, "wb");

    printf("plotting cqt spectrogram\n");

    int width = cqt_data[0].size();
    int height = cqt_data.size();

    // initialize row buffer
    png_bytep row = new png_byte[3 * width];
    PNGImage png_image(cqt_spectrogram_filepath, width, height, nullptr);

    double nmax = 0;
    double nmin = 0;
    double nrange = 0;
    for (auto &col : cqt_data) {
      for (double d : col) {
        nmax = std::max(d, nmax);
        nmin = std::min(d, nmin);
      }
    }
    nrange = nmax - nmin;

    for (auto &col : cqt_data) {
      fprintf(fdata, "%s", feature);

      for (int i = 0, j = 0; i < col.size(); ++i, ++j) {
        float level = (col[i] - nmin)/nrange;
        val2rgb2(level, row+(j*3) );

        fprintf(fdata, ",%f", level);
      }

      fprintf(fdata, "\n");

      png_image.WriteRow(row);
    }
    png_image.FinishWrite();

    delete [] row;
    fclose(fdata);

    printf("    image size: %dx%d\n", width, height);
  }

  fclose(file);

  printf("          input: %s\n", raw_pcm_file);
  printf("cqt spectrogram: %s\n", cqt_spectrogram_filepath);
  printf("       cqt data: %s\n", cqt_data_filepath);
}

int main(int argc, const char *argv[]) {
  //const char *raw_pcm_file = "/Users/neevek/Desktop/audioprocessing/note_scale/acou/23,24-acoustic-guitar.wav";
  //const char *raw_pcm_file = "/Users/neevek/Desktop/audioprocessing/note_scale/acou/30-acou.raw";

  //do_cqt(
      //"input/a.raw",
      //"out/a.png",
      //"out/a.txt",
      //"24"
      //);
  do_cqt(
      "/Users/neevek/Desktop/audioprocessing/input/Dataset_for_MIREX/WavfileRaw/abjones_5_09.raw",
      "/Users/neevek/Desktop/audioprocessing/input/Dataset_for_MIREX/WavfileCqtSpectrogram/abjones_5_09.png",
      "/Users/neevek/Desktop/audioprocessing/input/Dataset_for_MIREX/WavfileCqtData/abjones_5_09.txt",
      "24"
      );
  //do_cqt(
      //"input/28-acoustic-guitar.raw",
      //"out/28-acoustic-guitar.png",
      //"out/28-acoustic-guitar.txt",
      //"28"
      //);
  //do_cqt(
      //"input/29-acoustic-guitar.raw",
      //"out/29-acoustic-guitar.png",
      //"out/29-acoustic-guitar.txt",
      //"29"
      //);
  //do_cqt(
      //"input/31-acoustic-guitar.raw",
      //"out/31-acoustic-guitar.png",
      //"out/31-acoustic-guitar.txt",
      //"31"
      //);
  //do_cqt(
      //"input/33-acoustic-guitar.raw",
      //"out/33-acoustic-guitar.png",
      //"out/33-acoustic-guitar.txt",
      //"33"
      //);

  return 0;
}
