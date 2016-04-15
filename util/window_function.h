/**
 * No rights reserved, take it if it is useful to you.
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */
#ifndef WINDOW_FUNCTION_H_
#define WINDOW_FUNCTION_H_
#include <string>
#include <functional>

// http://en.wikipedia.org/wiki/Window_function
float hanning_window(int n, int N) {
  return 0.5f * (1-cos(2*M_PI*n/(N-1)));
}
float hamming_window(int n, int N) {
  return 0.54f - 0.46f * (float)cos((2 * M_PI * n) / (N - 1));
}
float nuttall_window(int n, int N) {
  return 0.355768 - 0.487396 * cos(2*M_PI*n/(N-1)) + 0.144232 * cos(4*M_PI*n/(N-1)) - 0.012604 * cos(6*M_PI*n/(N-1));
}
float blackman_window(int n, int N) {
  return 0.42659 - 0.49656 * cos(2*M_PI*n/(N-1)) + 0.076849 * cos(4*M_PI*n/(N-1));
}

float gaussian(float x, int N) {
  return pow(M_E, -pow((x-(N-1)/2)/(2*0.1f*N), 2));
}
float approximate_confined_gaussian_window(int n, int N) {
  return gaussian(n, N) - (gaussian(-0.5f, N)*(gaussian(n+N,N)+gaussian(n-N,N)))/
    (gaussian(-0.5f+N,N)+gaussian(-0.5f-N,N));
}

float gaussian_window(int n, int N) {
  return pow(M_E, -0.5f*pow((n-(N-1)/2)/(0.5f*(N-1)/2), 2));
}

float hann_poisson_window(int n, int N) {
  return 0.5f*(1-cos(2*M_PI*n)/(N-1))*pow(M_E, (-2*abs(N-1-2*n))/(N-1));
}

float rectangular_window(int n, int N) {
  return 1.f;
}

struct winfun_wrapper {
  const char *name;
  std::function<float(int,int)> winfun;
};

winfun_wrapper winfun_array[] = {
  {"rectangular", rectangular_window},
  {"hanning", hanning_window},
  {"hamming", hamming_window},
  {"nuttall", nuttall_window},
  {"blackman", blackman_window},
  {"approximate confined gaussian", approximate_confined_gaussian_window},
  {"gaussian", gaussian_window},
  {"hann poisson", hann_poisson_window},
};

std::string winfuns_to_string() {
  int size = sizeof(winfun_array)/sizeof(winfun_array[0]);

  std::string s;
  for (int i = 0; i < size; ++i) {
    s.append("\t   ");
    s.append(std::to_string(i));
    s.append(": ");
    s.append(winfun_array[i].name);
    s.append("\n");
  }

  return s;
}

#endif /* end of include guard: WINDOW_FUNCTION_H_ */
