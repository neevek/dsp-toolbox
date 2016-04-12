#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
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
  rgb[0] = rgb[1] = rgb[2] = 255 - fmin(level * 255, 255);
  //rgb[1] = fmin(level * 255, 255);
  //rgb[2] = fmin(level * 255, 255);
} 

// http://en.wikipedia.org/wiki/Window_function
float hamming_window(int n, int N) {
  return 0.54f - 0.46f * (float)cos((2 * M_PI * n) / (N - 1));
}

static void do_cqt(const char *raw_pcm_filepath,
                   const char *cqt_spectrogram_filepath,
                   const char *cqt_data_filepath,
                   unsigned int sample_rate,
                   float overlap,
                   double min_freq,
                   double max_freq,
                   unsigned int bins_per_octave,
                   bool is_stereo,
                   int use_channels,  // 0=both, 1=left, 2=right
                   const char *feature) {

  printf("\n");
  printf("================== debug info ==================\n");
  printf(" pcm input file: %s\n", raw_pcm_filepath);
  printf("cqt spectr file: %s\n", cqt_spectrogram_filepath);
  printf("  cqt data file: %s\n", cqt_data_filepath);
  printf("    sample rate: %d\n", sample_rate);
  printf("        overlap: %.4f\n", overlap);
  printf("  min frequency: %.4lf\n", min_freq);
  printf("  max frequency: %.4lf\n", max_freq);
  printf("bins per octave: %d\n", bins_per_octave);
  printf("      is stereo: %d\n", is_stereo);
  printf("   use channels: %d (0=both, 1=left, 2=right)\n", use_channels);
  printf("        feature: %s\n", feature);

  double sparse_cqt_kernel_threshold = 0.0054;

  Simac::ConstantQTransform cqt(sample_rate, min_freq, max_freq, bins_per_octave);
  unsigned int frame_size = cqt.getfftlength();
  FourierTransform fft(frame_size, 1, 0);
  cqt.sparsekernel(sparse_cqt_kernel_threshold);

  printf("        octaves: %d\n", cqt.getOctaves());
  printf("  est. min freq: %.4lf\n", cqt.getEstimatedMinFreq());

  FILE *file = fopen(raw_pcm_filepath, "rb");
  short inbuf[frame_size * 2];
  double frame_data[frame_size];

  std::vector<std::vector<double>> cqt_data;
  while(1) {
    if (is_stereo) {
      int n = fread(inbuf, sizeof(short) * 2, frame_size, file);
      if (n != frame_size) {
        break;
      }

      for (int i = 0; i < frame_size; ++i) {
        //frame_data[i] = inbuf[i * 2] + inbuf[i * 2 + 1];
        // use right channel only in this case
        //printf("i=%d, hamming: %f\n", i, hamming_window(i, frame_size));

        float hamming = hamming_window(i, frame_size);
        if (use_channels == 1) {
          frame_data[i] = inbuf[i * 2] * hamming;
        } else if (use_channels == 2) {
          frame_data[i] = inbuf[i * 2 + 1] * hamming;
        } else {
          frame_data[i] = (inbuf[i * 2] + inbuf[i * 2 + 1]) * hamming;
        }
      }

      long overlap_size = (long)(frame_size * sizeof(short) * 2 * overlap);
      if (overlap_size % 2 != 0) {
        --overlap_size;
      }
      fseek(file, -overlap_size, SEEK_CUR);
      //fseek(file, ftell(file)-overlap_size, SEEK_SET);
      //printf("overlap=%ld, frame_size=%d, cur_pos=%ld\n", overlap_size, frame_size, ftell(file));
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
  }

  if (cqt_data.size() > 0) {
    FILE *fdata = fopen(cqt_data_filepath, "wb");

    int width = cqt_data[0].size();
    int height = cqt_data.size();

    // initialize row buffer
    png_bytep row = new png_byte[3 * width];
    PNGImage png_image(cqt_spectrogram_filepath, width, height, nullptr);

    //double nmax = 0;
    //double nmin = 0;
    //double nrange = 0;
    //for (auto &col : cqt_data) {
      //for (double d : col) {
        //nmax = std::max(d, nmax);
        //nmin = std::min(d, nmin);
      //}
    //}
    //nrange = nmax - nmin;

    for (auto &col : cqt_data) {
      fprintf(fdata, "%s", feature);

      double nsum = 0;
      double nmax = 0;
      double nmin = 0;
      double nrange = 0;
      for (double d : col) {
        nmax = std::max(d, nmax);
        nmin = std::min(d, nmin);
        nsum += d;
      }
      nrange = nmax - nmin;

      for (int i = 0, j = 0; i < col.size(); ++i, ++j) {
        //float level = (col[i] - nmin)/nrange;
        //float level = pow(col[i] / nsum * 100, 2) / 100.f;
        float level = nsum == 0 ? 0 : col[i] / nsum;
        val2rgb2(level, row+(j*3) );

        fprintf(fdata, ",%f", level);
      }

      fprintf(fdata, "\n");

      png_image.WriteRow(row);
    }
    png_image.FinishWrite();

    delete [] row;
    fclose(fdata);

    printf("     image size: %dx%d\n", width, height);
  }

  fclose(file);

  printf("\n");
}

static char *raw_pcm_filepath = nullptr;
static char *cqt_spectrogram_filepath = nullptr;
static char *cqt_data_filepath = nullptr;
static bool is_stereo = false;
static int use_channels = 0; // 0=both, 1=left, 2=right
static unsigned int sample_rate = 0;
static float overlap = 0.5f;
static double min_freq = 246.94f; //61.735f;
static double max_freq = 4186.0f;
static unsigned int bins_per_octave = 36;
static char *feature = nullptr;

static void usage() {
  fprintf(stderr, "\nUsage:\n"
      "\t-i path to the raw PCM file\n"
      "\t-p path to output spectrogram image file\n"
      "\t-d path to cqt data text file\n"
      "\t-c is the raw PCM file stereo\n"
      "\t-u use which channels, default=%d (0=both, 1=left, 2=right)\n"
      "\t-s sample rate of the raw PCM file, default=%d\n"
      "\t-o overlap, default=%.4f\n"
      "\t-m min frequency, default=%.4lf\n"
      "\t-n max frequency, default=%.4lf\n"
      "\t-b bins per octave, default=%d\n"
      "\t-f feature flag\n"
      "\t-h print this help\n\n",
      use_channels,
      sample_rate,
      overlap,
      min_freq,
      max_freq,
      bins_per_octave
      );
}

static void log_and_exit(const char *msg) {
  fprintf(stderr, "%s", msg);
  usage();
  exit(1);
}

static void process_args(int argc, char **argv) {
  while (1) {
    int c = getopt(argc, argv, "i:p:d:c:u:s:o:m:n:b:f:h");
    if (c == -1) {
      break;
    }

    switch(c) {
      case 'i': raw_pcm_filepath = optarg; break;
      case 'p': cqt_spectrogram_filepath = optarg; break;
      case 'd': cqt_data_filepath = optarg; break;
      case 'c': is_stereo = (int)atoi(optarg); break;
      case 'u': use_channels = (int)atoi(optarg); break;
      case 's': sample_rate = (unsigned int)atoi(optarg); break;
      case 'o': overlap = (float )atof(optarg); break;
      case 'm': min_freq = (double)atof(optarg); break;
      case 'n': max_freq = (double)atof(optarg); break;
      case 'b': bins_per_octave = (int)atoi(optarg); break;
      case 'f': feature = optarg; break;
      case '?': usage(); break;
      default:
                usage();
                exit(1);
                break;
    }
  }
  if (raw_pcm_filepath == NULL) {
    log_and_exit("ERROR: must specify raw_pcm_filepath\n");
  }
  if (cqt_spectrogram_filepath == NULL) {
    log_and_exit("ERROR: must specify path to output spectrogram image file\n");
  }
  if (cqt_data_filepath == NULL) {
    log_and_exit("ERROR: must specify path to output cqt data file\n");
  }
  if (sample_rate <= 0) {
    log_and_exit("ERROR: must specify sample rate of the PCM file\n");
  }
  if (min_freq < 0) {
    log_and_exit("ERROR: must specify correct min_freq\n");
  }
  if (max_freq < 0) {
    log_and_exit("ERROR: must specify correct max_freq\n");
  }
  if (bins_per_octave <= 0 || (bins_per_octave % 12) != 0) {
    log_and_exit("ERROR: bins_per_octave must be multiples of 12\n");
  }
}

int main(int argc, const char *argv[]) {
  process_args(argc, (char **)argv);

  do_cqt(
      raw_pcm_filepath,
      cqt_spectrogram_filepath,
      cqt_data_filepath,
      sample_rate,
      overlap,
      min_freq,
      max_freq,
      bins_per_octave,
      is_stereo,
      use_channels,
      feature
      );

  return 0;
}
