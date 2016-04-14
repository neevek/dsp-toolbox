/**
 * No rights reserved, take it if it is useful to you.
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "ConstantQTransform.hxx"
#include "FourierTransform.hxx"
#include "png_image.h"
#include "window_function.h"
#include "util.h"

struct gen_cqt_cfg {
  std::string raw_pcm_filepath;
  std::string cqt_spectrogram_filepath;
  std::string cqt_data_filepath;
  bool is_stereo = false;
  int use_channels = 0; // 0=both, 1=left, 2=right
  unsigned int sample_rate = 0;
  float overlap = 0.5f;
  int winfun_type = 0; // default=hanning
  double min_freq = 123.47f; // 246.94f; //61.735f;
  double max_freq = 4186.0f;
  unsigned int bins_per_octave = 36;
  std::string feature;
  bool normalize_cqt_data = true;

  unsigned int frame_size;
};

static void write_spectrogram_and_cqt_data(
    const gen_cqt_cfg &cfg,
    const std::vector<std::vector<double>> &cqt_data) {
  if (cqt_data.size() == 0) {
    return;
  }

  int width = cqt_data[0].size();
  int height = cqt_data.size();

  PNGImage png_image(cfg.cqt_spectrogram_filepath, width, height, nullptr);

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

  FILE *fdata = fopen(cfg.cqt_data_filepath.c_str(), "wb");
  fprintf(fdata, "%d,%d,%.4f\n", cfg.sample_rate, cfg.frame_size, cfg.overlap);
  // initialize row buffer
  png_bytep row = new png_byte[3 * width];

  for (auto &col : cqt_data) {
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

      fprintf(fdata, ",%f", cfg.normalize_cqt_data ? level : col[i]);
    }

    fprintf(fdata, "\n");

    png_image.WriteRow(row);
  }
  delete [] row;
  png_image.Finish();

  fclose(fdata);

  fprintf(stderr, "     image size: %dx%d\n", width, height);
}

static void do_cqt(gen_cqt_cfg cfg) {
  double sparse_cqt_kernel_threshold = 0.0054;

  Simac::ConstantQTransform cqt(cfg.sample_rate, cfg.min_freq, cfg.max_freq, cfg.bins_per_octave);
  cfg.frame_size = cqt.getfftlength();
  FourierTransform fft(cfg.frame_size, 1, 0);
  cqt.sparsekernel(sparse_cqt_kernel_threshold);

  winfun_wrapper &selected_winfun = winfun_array[cfg.winfun_type];

  fprintf(stderr, "\n");
  fprintf(stderr, "================== debug info ==================\n");
  fprintf(stderr, " pcm input file: %s\n", cfg.raw_pcm_filepath.c_str());
  fprintf(stderr, "cqt spectr file: %s\n", cfg.cqt_spectrogram_filepath.c_str());
  fprintf(stderr, "  cqt data file: %s\n", cfg.cqt_data_filepath.c_str());
  fprintf(stderr, "    sample rate: %d\n", cfg.sample_rate);
  fprintf(stderr, "        overlap: %.4f\n", cfg.overlap);
  fprintf(stderr, "window function: %s\n", selected_winfun.name);
  fprintf(stderr, "  min frequency: %.4lf\n", cfg.min_freq);
  fprintf(stderr, "  max frequency: %.4lf\n", cfg.max_freq);
  fprintf(stderr, "bins per octave: %d\n", cfg.bins_per_octave);
  fprintf(stderr, "      is stereo: %d\n", cfg.is_stereo);
  fprintf(stderr, "   use channels: %d (0=both, 1=left, 2=right)\n", cfg.use_channels);
  fprintf(stderr, "        feature: %s\n", cfg.feature.c_str());
  fprintf(stderr, "        octaves: %d\n", cqt.getOctaves());
  fprintf(stderr, "  est. min freq: %.4lf\n", cqt.getEstimatedMinFreq());
  fprintf(stderr, "     frame size: %d\n", cfg.frame_size);
  fprintf(stderr," cqt normalized: %d\n", cfg.normalize_cqt_data);

  FILE *file = fopen(cfg.raw_pcm_filepath.c_str(), "rb");
  short inbuf[cfg.frame_size * 2];
  double frame_data[cfg.frame_size];

  std::vector<std::vector<double>> cqt_data;
  while(1) {
    if (cfg.is_stereo) {
      int n = fread(inbuf, sizeof(short) * 2, cfg.frame_size, file);
      if (n != cfg.frame_size) {
        break;
      }

      for (int i = 0; i < cfg.frame_size; ++i) {
        float window = selected_winfun.winfun(i, cfg.frame_size);
        if (cfg.use_channels == 1) {
          frame_data[i] = inbuf[i * 2] * window;
        } else if (cfg.use_channels == 2) {
          frame_data[i] = inbuf[i * 2 + 1] * window;
        } else {
          frame_data[i] = (inbuf[i * 2] + inbuf[i * 2 + 1]) * window;
        }
      }

      long overlap_size = (long)(cfg.frame_size * sizeof(short) * 2 * cfg.overlap);
      if (overlap_size % 2 != 0) {
        --overlap_size;
      }
      fseek(file, -overlap_size, SEEK_CUR);
    } else {
      int n = fread(inbuf, sizeof(short), cfg.frame_size, file);
      if (n != cfg.frame_size) {
        break;
      }

      for (int i = 0; i < cfg.frame_size; ++i) {
        float window = selected_winfun.winfun(i, cfg.frame_size);
        frame_data[i] = inbuf[i] * window;
      }

      long overlap_size = (long)(cfg.frame_size * sizeof(short) * cfg.overlap);
      if (overlap_size % 2 != 0) {
        --overlap_size;
      }
      fseek(file, -overlap_size, SEEK_CUR);
    }

    fft.doIt(frame_data);
    cqt.doIt(fft.spectrum());

    std::vector<double> spectrum = cqt.constantQSpectrum();
    std::vector<double> col(spectrum);

    for (int i = 0; i < col.size()/2; ++i) {
      col[i] = sqrt(col[i*2] * col[i*2] + col[i*2+1] * col[i*2+1]);
    }
    col.resize(col.size()/2);

    cqt_data.push_back(std::move(col));
  }
  fclose(file);

  write_spectrogram_and_cqt_data(cfg, cqt_data);

  fprintf(stderr, "\n");
}

static void usage(const gen_cqt_cfg &cfg) {
  fprintf(stderr, "\nUsage:\n"
      "\t-i path to the raw PCM file\n"
      "\t-p path to output spectrogram image file\n"
      "\t-d path to cqt data text file\n"
      "\t-c is the raw PCM file stereo\n"
      "\t-u use which channels, default=%d (0=both, 1=left, 2=right)\n"
      "\t-s sample rate of the raw PCM file, default=%d\n"
      "\t-o overlap, default=%.4f\n"
      "\t-w window function, one of the following (default=0): \n"
      "%s"
      "\t-m min frequency, default=%.4lf\n"
      "\t-n max frequency, default=%.4lf\n"
      "\t-b bins per octave, default=%d\n"
      "\t-f feature flag\n"
      "\t-x normalize cqt data (default=1)\n"
      "\t-h print this help\n\n",
      cfg.use_channels,
      cfg.sample_rate,
      cfg.overlap,
      winfuns_to_string().c_str(),
      cfg.min_freq,
      cfg.max_freq,
      cfg.bins_per_octave
      );
}

static void log_and_exit(const gen_cqt_cfg &cfg, const char *msg) {
  fprintf(stderr, "%s", msg);
  usage(cfg);
  exit(1);
}

static gen_cqt_cfg process_args(int argc, char **argv) {
  gen_cqt_cfg cfg;

  while (1) {
    int c = getopt(argc, argv, "i:p:d:c:u:s:o:w:m:n:b:f:x:h");
    if (c == -1) {
      break;
    }

    switch(c) {
      case 'i': cfg.raw_pcm_filepath = optarg; break;
      case 'p': cfg.cqt_spectrogram_filepath = optarg; break;
      case 'd': cfg.cqt_data_filepath = optarg; break;
      case 'c': cfg.is_stereo = (int)atoi(optarg); break;
      case 'u': cfg.use_channels = (int)atoi(optarg); break;
      case 's': cfg.sample_rate = (unsigned int)atoi(optarg); break;
      case 'o': cfg.overlap = (float )atof(optarg); break;
      case 'w': cfg.winfun_type = atoi(optarg); break;
      case 'm': cfg.min_freq = (double)atof(optarg); break;
      case 'n': cfg.max_freq = (double)atof(optarg); break;
      case 'b': cfg.bins_per_octave = (int)atoi(optarg); break;
      case 'f': cfg.feature = optarg; break;
      case 'x': cfg.normalize_cqt_data = (int)atoi(optarg); break;
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
  if (cfg.cqt_spectrogram_filepath.size() == 0) {
    log_and_exit(cfg, "ERROR: must specify path to output spectrogram image file\n");
  }
  if (cfg.cqt_data_filepath.size() == 0) {
    log_and_exit(cfg, "ERROR: must specify path to output cqt data file\n");
  }
  if (cfg.sample_rate <= 0) {
    log_and_exit(cfg, "ERROR: must specify sample rate of the PCM file\n");
  }
  if (cfg.min_freq < 0) {
    log_and_exit(cfg, "ERROR: must specify correct min_freq\n");
  }
  if (cfg.max_freq < 0) {
    log_and_exit(cfg, "ERROR: must specify correct max_freq\n");
  }
  if (cfg.bins_per_octave <= 0 || (cfg.bins_per_octave % 12) != 0) {
    log_and_exit(cfg, "ERROR: bins_per_octave must be multiples of 12\n");
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
  do_cqt(process_args(argc, (char **)argv));
  return 0;
}
