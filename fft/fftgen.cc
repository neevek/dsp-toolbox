/**
 * No rights reserved, take it if it is useful to you.
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */
#include <cstdio>
#include <unistd.h>
#include <cstdlib>
#include <vector>
#include "kiss_fft.h"
#include "kiss_fftr.h"
#include "../util/png_image.h"
#include "../util/window_function.h"
#include "../util/util.h"

struct gen_fft_cfg {
  std::string raw_pcm_filepath;
  std::string fft_spectrogram_filepath;
  std::string fft_data_filepath;
  bool is_stereo = false;
  int use_channels = 0; // 0=both, 1=left, 2=right
  unsigned int sample_rate = 0;
  int herz_per_bin = 5;
  int reserved_freq_bin = 200;
  float overlap = 0.5f;
  int winfun_type = 0; // default=hanning
  std::string feature;
  bool normalize_fft_data = true;
  bool add_redundancy_zero = false;

  unsigned int frame_size;
};

static void write_spectrogram_and_fft_data(
    const gen_fft_cfg &cfg,
    const std::vector<std::vector<float>> &fft_data) {
  if (fft_data.size() == 0) {
    return;
  }

  int width = fft_data[0].size();
  int height = fft_data.size();

  PNGImage png_image(cfg.fft_spectrogram_filepath, width, height, nullptr);

  //double nmax = 0;
  //double nmin = 0;
  //double nrange = 0;
  //for (auto &col : fft_data) {
  //for (double d : col) {
  //nmax = std::max(d, nmax);
  //nmin = std::min(d, nmin);
  //}
  //}
  //nrange = nmax - nmin;

  FILE *fdata = fopen(cfg.fft_data_filepath.c_str(), "wb");
  fprintf(fdata, "%d,%d,%.4f\n", cfg.sample_rate, cfg.frame_size, cfg.overlap);
  // initialize row buffer
  png_bytep row = new png_byte[3 * width];

  for (auto &col : fft_data) {
    fprintf(fdata, "%s", cfg.feature.c_str());

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
      float level = nsum == 0 ? 0 : col[i] / nsum;
      val_to_rgb(level, row+(j*3) );

      fprintf(fdata, ",%f", cfg.normalize_fft_data ? level : col[i]);
    }

    fprintf(fdata, "\n");

    png_image.WriteRow(row);
  }
  delete [] row;
  png_image.Finish();

  fclose(fdata);

  fprintf(stderr, "     image size: %dx%d\n", width, height);
}

static void do_fft(gen_fft_cfg cfg) {
  int nfft = cfg.sample_rate / cfg.herz_per_bin;
  int nfreqs = nfft/2 + 1;
  int restricted_nfreqs = fmin(nfreqs, cfg.reserved_freq_bin);

  cfg.frame_size = nfft;

  kiss_fftr_cfg fft_cfg = kiss_fftr_alloc(nfft, 0, 0, 0);
  kiss_fft_scalar *fft_input = (kiss_fft_scalar *)malloc(sizeof(kiss_fft_scalar) * nfft);
  kiss_fft_cpx *fft_output = (kiss_fft_cpx *)malloc(sizeof(kiss_fft_cpx) * nfreqs);
  short *inbuf = (short *)malloc(sizeof(short) * nfft * 2); // for stereo

  winfun_wrapper &selected_winfun = winfun_array[cfg.winfun_type];

  fprintf(stderr, "\n");
  fprintf(stderr, "================== debug info ==================\n");
  fprintf(stderr, " pcm input file: %s\n", cfg.raw_pcm_filepath.c_str());
  fprintf(stderr, "fft spectr file: %s\n", cfg.fft_spectrogram_filepath.c_str());
  fprintf(stderr, "  fft data file: %s\n", cfg.fft_data_filepath.c_str());
  fprintf(stderr, "    sample rate: %d\n", cfg.sample_rate);
  fprintf(stderr, "        overlap: %.4f\n", cfg.overlap);
  fprintf(stderr, "window function: %s\n", selected_winfun.name);
  fprintf(stderr, "      is stereo: %d\n", cfg.is_stereo);
  fprintf(stderr, "   use channels: %d (0=both, 1=left, 2=right)\n", cfg.use_channels);
  fprintf(stderr, "        feature: %s\n", cfg.feature.c_str());
  fprintf(stderr, "     Hz per bin: %d\n", cfg.herz_per_bin);
  fprintf(stderr, "  reserved bins: %d\n", cfg.reserved_freq_bin);
  fprintf(stderr, "     frame size: %d\n", cfg.frame_size);
  fprintf(stderr, " fft normalized: %d\n", cfg.normalize_fft_data);
  fprintf(stderr, "       add zero: %d\n", cfg.add_redundancy_zero);

  FILE *file = fopen(cfg.raw_pcm_filepath.c_str(), "rb");

  int arz = cfg.add_redundancy_zero;

  std::vector<std::vector<float>> fft_data;
  while (1) {
    if (cfg.is_stereo) {
      int n = fread(inbuf, sizeof(short) * 2, cfg.frame_size / (arz ? 2 : 1), file);
      if ((!arz && n != cfg.frame_size) || (arz && n != cfg.frame_size/2)) {
        break;
      }

      for (int i = 0; i < cfg.frame_size; ++i) {
        float window = selected_winfun.winfun(i, cfg.frame_size);
        if (cfg.use_channels == 1) {
          fft_input[i] = inbuf[i * 2] * window;
        } else if (cfg.use_channels == 2) {
          fft_input[i] = inbuf[i * 2 + 1] * window;
        } else {
          fft_input[i] = (inbuf[i * 2] + inbuf[i * 2 + 1]) * window;
        }
      }

      long overlap_size = (long)(cfg.frame_size * sizeof(short) * 2 * cfg.overlap);
      if (overlap_size % 2 != 0) {
        --overlap_size;
      }
      fseek(file, -overlap_size, SEEK_CUR);
    } else {
      int n = fread(inbuf, sizeof(short), cfg.frame_size / (arz ? 2 : 1), file);
      if ((!arz && n != cfg.frame_size) || (arz && n != cfg.frame_size/2)) {
        break;
      }

      for (int i = 0; i < cfg.frame_size; ++i) {
        float window = selected_winfun.winfun(i, cfg.frame_size);
        fft_input[i] = inbuf[i] * window;
      }

      long overlap_size = (long)(cfg.frame_size * sizeof(short) * cfg.overlap);
      if (overlap_size % 2 != 0) {
        --overlap_size;
      }
      fseek(file, -overlap_size, SEEK_CUR);
    }

    if (arz) {
      for (int i = cfg.frame_size/2-1, j = cfg.frame_size-1; i > 0; --i, j-=2) {
        fft_input[j] = 0;
        fft_input[j-1] = fft_input[i];
      }
    }

    /*if (remove_dc) {*/
    /*float avg = 0;*/
    /*for (int i=0;i<nfft;++i)  avg += fft_input[i];*/
    /*avg /= nfft;*/
    /*for (int i=0;i<nfft;++i)  fft_input[i] -= (kiss_fft_scalar)avg;*/
    /*}*/

    kiss_fftr(fft_cfg, fft_input, fft_output);

    std::vector<float> col;

    for (int i = 0; i < restricted_nfreqs; ++i) {
      col.push_back(sqrt(fft_output[i].r * fft_output[i].r + fft_output[i].i * fft_output[i].i));
    }

    fft_data.push_back(std::move(col));
  }

  free(fft_cfg);
  free(fft_input);
  free(fft_output);
  free(inbuf);
  fclose(file);

  write_spectrogram_and_fft_data(cfg, fft_data);

  fprintf(stderr, "\n");
}

static void usage(const gen_fft_cfg &cfg) {
  fprintf(stderr, "\nUsage:\n"
      "\t-i path to the raw PCM file\n"
      "\t-p path to output spectrogram image file\n"
      "\t-d path to fft data text file\n"
      "\t-c is the raw PCM file stereo\n"
      "\t-u use which channels, default=%d (0=both, 1=left, 2=right)\n"
      "\t-s sample rate of the raw PCM file, default=%d\n"
      "\t-o overlap, default=%.4f\n"
      "\t-w window function, one of the following (default=0): \n"
      "%s"
      "\t-m Hz per bin, default=%d\n"
      "\t-n reserved bins, default=%d\n"
      "\t-f feature flag\n"
      "\t-x normalize fft data (default=1)\n"
      "\t-z add redundancy zero (default=0)\n"
      "\t-h print this help\n\n",
      cfg.use_channels,
      cfg.sample_rate,
      cfg.overlap,
      winfuns_to_string().c_str(),
      cfg.herz_per_bin,
      cfg.reserved_freq_bin
      );
}

static void log_and_exit(const gen_fft_cfg &cfg, const char *msg) {
  fprintf(stderr, "%s", msg);
  usage(cfg);
  exit(1);
}

static gen_fft_cfg process_args(int argc, char **argv) {
  gen_fft_cfg cfg;

  while (1) {
    int c = getopt(argc, argv, "i:p:d:c:u:s:o:w:m:n:f:x:z:h");
    if (c == -1) {
      break;
    }

    switch(c) {
      case 'i': cfg.raw_pcm_filepath = optarg; break;
      case 'p': cfg.fft_spectrogram_filepath = optarg; break;
      case 'd': cfg.fft_data_filepath = optarg; break;
      case 'c': cfg.is_stereo = (int)atoi(optarg); break;
      case 'u': cfg.use_channels = (int)atoi(optarg); break;
      case 's': cfg.sample_rate = (unsigned int)atoi(optarg); break;
      case 'o': cfg.overlap = (float )atof(optarg); break;
      case 'w': cfg.winfun_type = atoi(optarg); break;
      case 'm': cfg.herz_per_bin= (double)atof(optarg); break;
      case 'n': cfg.reserved_freq_bin= (double)atof(optarg); break;
      case 'f': cfg.feature = optarg; break;
      case 'x': cfg.normalize_fft_data = (int)atoi(optarg); break;
      case 'z': cfg.add_redundancy_zero = (int)atoi(optarg); break;
      case '?': usage(cfg); break;
      default:
                usage(cfg);
                exit(1);
                break;
    }
  }
  if (cfg.raw_pcm_filepath.size() == 0) {
    log_and_exit(cfg, "ERROR: must specify raw_pcm_filepath\n");
  }
  if (cfg.fft_spectrogram_filepath.size() == 0) {
    log_and_exit(cfg, "ERROR: must specify path to output spectrogram image file\n");
  }
  if (cfg.fft_data_filepath.size() == 0) {
    log_and_exit(cfg, "ERROR: must specify path to output fft data file\n");
  }
  if (cfg.sample_rate <= 0) {
    log_and_exit(cfg, "ERROR: must specify sample rate of the PCM file\n");
  }
  if (cfg.herz_per_bin <= 0) {
    log_and_exit(cfg, "ERROR: must specify correct herz_per_bin\n");
  }
  if (cfg.reserved_freq_bin <= 0) {
    log_and_exit(cfg, "ERROR: must specify correct reserved_freq_bin\n");
  }

  int size = sizeof(winfun_array)/sizeof(winfun_array[0]);
  if (cfg.winfun_type < 0 || cfg.winfun_type >= size) {
    log_and_exit(cfg, "ERROR: window function not found\n");
  }

  if (cfg.winfun_type == 0) {
    cfg.overlap = 0;  // if no window function is selected, overlap will be ignored
  }

  return cfg;
}

int main(int argc, const char *argv[]) {
  do_fft(process_args(argc, (char **)argv));
  return 0;
}
