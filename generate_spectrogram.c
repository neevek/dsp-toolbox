#include <stdio.h>
#include <math.h>
#include <png.h>
#include <unistd.h>

#include "kiss_fft.h"
#include "kiss_fftr.h"

float minval = 0;
float maxval = 0;
float valrange = 0;

#define CHECKNULL(p) if ( (p)==NULL ) do { fprintf(stderr,"CHECKNULL failed @ %s(%d): %s\n",__FILE__,__LINE__,#p );exit(1);} while(0)

/*static void val2rgb(float x, png_byte *p) {*/
    /*const double pi = 3.14159265358979;*/
    /*p[0] = (int)(255*sin(x*pi));*/
    /*p[1] = (int)(255*fabs(sin(x*pi*3/2)));*/
    /*p[2] = (int)(255*fabs(sin(x*pi*5/2)));*/
    /*//fprintf(stderr,"%.2f : %d,%d,%d\n",x,(int)p->r,(int)p->g,(int)p->b);*/
/*}*/

// Modified version of Dan Bruton's algorithm:
// http://www.physics.sfasu.edu/astro/color/spectra.html
static void val2rgb2(float level, png_byte *rgb) {
  rgb[0] = level * 255;
  rgb[1] = level * 255;
  rgb[2] = level * 255;
}

static void val2rgb3(float level, png_byte *rgb) {
  printf("%f\n", level);
    level *= 0.6625;
    double r = 0.0, g = 0.0, b = 0.0;
    if (level >= 0 && level < 0.15) {
        r = (0.15 - level) / (0.15 + 0.075);
        g = 0.0;
        b = 1.0;
    } else if (level >= 0.15 && level < 0.275) {
        r = 0.0;
        g = (level - 0.15) / (0.275 - 0.15);
        b = 1.0;
    } else if (level >= 0.275 && level < 0.325) {
        r = 0.0;
        g = 1.0;
        b = (0.325 - level) / (0.325 - 0.275);
    } else if (level >= 0.325 && level < 0.5) {
        r = (level - 0.325) / (0.5 - 0.325);
        g = 1.0;
        b = 0.0;
    } else if (level >= 0.5 && level < 0.6625) {
        r = 1.0;
        g = (0.6625 - level) / (0.6625 - 0.5f);
        b = 0.0;
    }

    // Intensity correction.
    double cf = 1.0;
    if (level >= 0.0 && level < 0.1) {
        cf = level / 0.1;
    }
    cf *= 255.0;

    // Pack RGB values into a 32-bit uint.
    rgb[0] = (r * cf + 0.5);
    rgb[1] = (g * cf + 0.5);
    rgb[2] = (b * cf + 0.5);
}

static void cpx2pixels(png_byte *res, const float *fbuf, size_t n) {
    size_t i;
    float minval,maxval,valrange;
    minval=maxval=fbuf[0];

    for (i = 0; i < n; ++i) {
        if (fbuf[i] > maxval) maxval = fbuf[i];
        if (fbuf[i] < minval) minval = fbuf[i];
    }

    fprintf(stderr,"min ==%f,max=%f\n",minval,maxval);
    valrange = maxval-minval;
    if (valrange == 0) {
        fprintf(stderr,"min == max == %f\n",minval);
        exit (1);
    }

    for (i = 0; i < n; ++i)
        val2rgb2( (fbuf[i] - minval)/valrange , res+(i*3) );
}

// not efficient enough, but this is not the bottle neck
static int is_sound_quiet(float *magbuf, int count, float threshold) {
  float valrange = maxval - minval;
  double total = 0;
  for (int i = 0; i < count; ++i) {
    total += magbuf[i];
  }
  return ((total / count) - minval)/valrange <= threshold;
}


int write_image(const char *filename, int width, int height, int nrows, float *buffer, char* title, float mag_to_ignore_threshold) {
  int code = 0;
  FILE *fp = NULL;
  png_structp png_ptr = NULL;
  png_infop info_ptr = NULL;
  png_bytep row = NULL;
  
  // Open file for writing (binary mode)
  fp = fopen(filename, "wb");
  if (fp == NULL) {
    fprintf(stderr, "Could not open file %s for writing\n", filename);
    code = 1;
    goto finalise;
  }

  // Initialize write structure
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL) {
    fprintf(stderr, "Could not allocate write struct\n");
    code = 1;
    goto finalise;
  }

  // Initialize info structure
  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL) {
    fprintf(stderr, "Could not allocate info struct\n");
    code = 1;
    goto finalise;
  }

  // Setup Exception handling
  if (setjmp(png_jmpbuf(png_ptr))) {
    fprintf(stderr, "Error during png creation\n");
    code = 1;
    goto finalise;
  }

  png_init_io(png_ptr, fp);

  // Write header (8 bit colour depth)
  png_set_IHDR(png_ptr, info_ptr, width, height,
      8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
      PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

  // Set title
  if (title != NULL) {
    png_text title_text;
    title_text.compression = PNG_TEXT_COMPRESSION_NONE;
    title_text.key = "Title";
    title_text.text = title;
    png_set_text(png_ptr, info_ptr, &title_text, 1);
  }

  png_write_info(png_ptr, info_ptr);

  // Allocate memory for one row (3 bytes per pixel - RGB)
  row = (png_bytep) malloc(3 * width * sizeof(png_byte));
  float valrange = maxval - minval;

  // Write image data
  for (int i = 0; i < nrows; ++i) {
    if (is_sound_quiet(&buffer[i * width], width, mag_to_ignore_threshold)) {
      continue;
    }

    for (int j = 0; j < width; ++j) {
      val2rgb2((buffer[i * width + j] - minval)/valrange , row+(j*3) );
    }
    png_write_row(png_ptr, row);
  }

  // End write
  png_write_end(png_ptr, NULL);

  finalise:
  if (fp != NULL) fclose(fp);
  if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
  if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
  if (row != NULL) free(row);

  return code;
}

static int write_magnitudes_to_file(
    const char *outfile, 
    float *magbuf,
    int nfreqs, 
    int nrows,
    int feature_flag,
    float mag_to_ignore_threshold
    ) {

  FILE *file = fopen(outfile, "wb");

  float valrange = maxval - minval;

  int lines_ignored = 0;
  for (int i = 0; i < nrows; ++i) {
    if (is_sound_quiet(&magbuf[i * nfreqs], nfreqs, mag_to_ignore_threshold)) {
      ++lines_ignored;
      continue;
    }

    fprintf(file, "%d,", feature_flag);

    for (int j = 0; j < nfreqs; ++j) {
      fprintf(file, "%f", (magbuf[i * nfreqs + j] - minval)/valrange);
      if (j + 1 == nfreqs) {
        fprintf(file, "\n");
      } else {
        fprintf(file, ",");
      }
    }
  }

  fclose(file);

  return lines_ignored;
}

static void do_fft(
    const char *audio_file_path, 
    const char *out_image_file,
    const char *out_magnitude_file,
    int herz_per_bin, 
    int reserved_freq_bin,
    int sample_rate,
    int stereo,
    int round_to_avg,
    int feature_flag // 1 or 0, being printed to the FFT data file
    ) {

  int nfft = sample_rate / herz_per_bin;
  int nfreqs = nfft/2 + 1;
  int restricted_nfreqs = fmin(nfreqs, reserved_freq_bin);

  kiss_fftr_cfg cfg = kiss_fftr_alloc(nfft, 0, 0, 0);
  kiss_fft_scalar *fft_input = (kiss_fft_scalar *)malloc(sizeof(kiss_fft_scalar) * nfft);
  kiss_fft_cpx *fft_output = (kiss_fft_cpx *)malloc(sizeof(kiss_fft_cpx) * nfreqs);
  short *inbuf = (short *)malloc(sizeof(short) * nfft * 2); // for stereo
  float *magbuf = (float *)malloc(sizeof(float) * restricted_nfreqs);
  memset(magbuf, 0, sizeof(magbuf[0]) * restricted_nfreqs);

  FILE *file = fopen(audio_file_path, "rb");

  int round_ctr = 0;

  float *img_data = NULL;
  int nrows = 0;

  while (1) {
    if (stereo) {
      int n = fread(inbuf, sizeof(short) * 2, nfft, file);
      if (n != nfft) {
        break;
      }

      for (int i = 0; i < nfft; ++i) {
        fft_input[i] = inbuf[i * 2] + inbuf[i * 2 + 1];
      }
    } else {
      int n = fread(inbuf, sizeof(short), nfft, file);
      if (n != nfft) {
        break;
      }

      for (int i = 0; i < nfft; ++i) {
        fft_input[i] = inbuf[i];
      }
    }

    /*if (remove_dc) {*/
    /*float avg = 0;*/
    /*for (int i=0;i<nfft;++i)  avg += fft_input[i];*/
    /*avg /= nfft;*/
    /*for (int i=0;i<nfft;++i)  fft_input[i] -= (kiss_fft_scalar)avg;*/
    /*}*/

    kiss_fftr(cfg, fft_input, fft_output);

    for (int i = 0; i < restricted_nfreqs; ++i) {
      magbuf[i] += fft_output[i].r * fft_output[i].r + fft_output[i].i * fft_output[i].i;
    }

    if (++round_ctr == round_to_avg) {
      round_ctr = 0;
      ++nrows;

      img_data = (float*)realloc(img_data, sizeof(float)*nrows*restricted_nfreqs);
      const float eps = 1;
      for (int i = 0; i < restricted_nfreqs; ++i) {
        /*img_data[(nrows - 1) * nfreqs + i] = magbuf[i];*/
        int index = (nrows - 1) * restricted_nfreqs + i;
        img_data[index] = 10 * log10(magbuf[i] / round_to_avg + eps);

        if (img_data[index] > maxval) {
          maxval = img_data[index];
        }
        if (img_data[index] < minval) {
          minval = img_data[index];
        }
      }
      memset(magbuf, 0, sizeof(magbuf[0]) * restricted_nfreqs);
    }
  }

  float mag_to_ignore_threshold = 0.0f;
  int lines_ignored = write_magnitudes_to_file(out_magnitude_file, img_data, restricted_nfreqs, nrows, feature_flag, mag_to_ignore_threshold);
  printf("%d, %d\n", lines_ignored, nrows);
  write_image(out_image_file, restricted_nfreqs, nrows - lines_ignored, nrows, img_data, NULL, mag_to_ignore_threshold);

  free(fft_input);
  free(fft_output);
  free(inbuf);
  free(magbuf);
  free(img_data);
  fclose(file);
}

static void usage() {
  fprintf(stderr, "Usage:\n"
      "\t-i s: path to the raw PCM file\n"
      "\t-p s: path to output spectrogram image file\n"
      "\t-d s: path to ttf data text file\n"
      "\t-z d: number of Herz per bin in frequency domain, default: 5\n"
      "\t-b d: number of bins to reserve in frequency domain, default: 200\n"
      "\t-s d: sample rate of the raw PCM file\n"
      "\t-c d: is the raw PCM file stereo\n"
      "\t-f d: feature flag\n"
      "\t-h: print this help info\n"
      );
}

static void log_and_exit(const char *msg) {
  fprintf(stderr, "%s", msg);
  exit(1);
}

static char *pcm_file = NULL;
static char *specgrm_file = NULL;
static char *fft_data_file = NULL;
static int herz_per_bin = 5;
static int bins_to_reserve = 200;
static int sample_rate = 0;  // must greater than 0
static int is_stereo = -1;   
static int round_to_avg = 1; 
static int feature_flag = -1; 

static void process_args(int argc, char **argv) {
  while (1) {
    int c = getopt(argc, argv, "i:p:d:z:b:s:c:r:f:h");
    if (c == -1) {
      break;
    }

    switch(c) {
      case 'i': pcm_file = optarg; break;
      case 'p': specgrm_file = optarg; break;
      case 'd': fft_data_file = optarg; break;
      case 'z': herz_per_bin = (int)atoi(optarg); break;
      case 'b': bins_to_reserve = (int)atoi(optarg); break;
      case 's': sample_rate = (int)atoi(optarg); break;
      case 'c': is_stereo = (int)atoi(optarg); break;
      case 'r': round_to_avg = (int)atoi(optarg); break;
      case 'f': feature_flag = (int)atoi(optarg); break;
      case '?': usage(); break;
      default:
                usage();
                exit(1);
                break;
    }
  }
  if (pcm_file == NULL) {
    log_and_exit("ERROR: must specify pcm_file\n");
  }
  if (specgrm_file == NULL) {
    log_and_exit("ERROR: must specify path to output spectrogram image file\n");
  }
  if (fft_data_file == NULL) {
    log_and_exit("ERROR: must specify path to output fft data file\n");
  }
  if (sample_rate == 0) {
    log_and_exit("ERROR: must specify sample rate of the PCM file\n");
  }
  if (is_stereo == -1) {
    log_and_exit("ERROR: must specify number of is_stereo of the PCM file\n");
  }
  if (feature_flag == -1) {
    log_and_exit("ERROR: must specify number of is_stereo of the PCM file\n");
  }
}

int main(int argc, const char *argv[]) {
  process_args(argc, (char **)argv);
  do_fft(pcm_file, specgrm_file, fft_data_file, herz_per_bin, bins_to_reserve, sample_rate, is_stereo, round_to_avg, feature_flag);

  /*do_fft("/Users/neevek/Desktop/guitarsound/do-2-cut.raw", 5, 200, "out/1.png", "out/1.txt");*/

  return 0;
}
