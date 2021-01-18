// Assignment 2
//
// XXX Analysis!
// - find the most prominent N peaks in a sequence of sound clips
// - break the input into a sequence of overlaping sound clips
//   + each sound clip is 2048 samples long
//   + each sound clip should overlap the last by 1024 samples (50% overlap)
//   + use a Hann window to prepare the data before FFT analysis
// - use the FFT to analyse each sound clip
//   + use an FFT size of 8192 (pad the empty data with zeros)
//   + find the frequency and amplitude of each peak in the clip
//     * find maxima, calculate frequency, sort these by amplitude
//   + store the frequency and amplitude of the most prominent N of these
// - test this whole process using several given audio files (e.g.,
// sine-sweep.wav, impulse-sweep.wav, sawtooth-sweep.wav)
//
// -- Karl Yerkes / 2021-01-12 / MAT240B
//

#include <algorithm>  // std::sort
#include <vector>
#include <complex>
#include <iostream>
#include <valarray>

const double PI = 3.141592653589793238460;
const int WINDOW_LENGTH = 2048;
const int PADDED_LENGTH = 8192;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
typedef std::pair<double,double> freq_amp;

// Cooley–Tukey FFT (in-place, divide-and-conquer)
// Higher memory requirements and redundancy although more intuitive
void fft(CArray& x) {
  const size_t N = x.size();
  //std::cout << "N = " << N << std::endl;
  if (N <= 1) return;

  // divide
  CArray even = x[std::slice(0, N / 2, 2)];
  CArray odd = x[std::slice(1, N / 2, 2)];

  // conquer
  //std::cout << "even" << std::endl;
  fft(even);
  //std::cout << "odd" << std::endl;
  fft(odd);

  // combine
  //std::cout << "combine" << std::endl;
  for (size_t k = 0; k < N / 2; ++k) {
    Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
    x[k] = even[k] + t;
    x[k + N / 2] = even[k] - t;
  }
  
  // std::cout << " N = " << N << " x.size() = " << x.size() << std::endl;
  // for (int i = 0; i < x.size(); i++) {
  //   std::cout << x[i] << std::endl;
  // }
}

// inverse fft (in-place)
void ifft(CArray& x) {
  // conjugate the complex numbers
  x = x.apply(std::conj);

  // forward fft
  fft(x);

  // conjugate the complex numbers again
  x = x.apply(std::conj);

  // scale the numbers
  x /= x.size();
}

void print_data(CArray& d) { // output vals
  std::cout << "size = " << d.size() << std::endl;
  for (int i = 0; i < d.size(); ++i) {
    //printf("%lf\n", d[i]);
    std::cout << d[i] << std::endl;
  }
}

double hann(double x) { // Hanning window
  return (1 - cos(2 * PI * (x - 1) / 2)) / 2;
}

// bool compare (Complex i,Complex j) { return (i.real()<j.real()); } // from http://www.cplusplus.com/reference/algorithm/sort/

int main(int argc, char *argv[]) {
  // take the data in
  std::vector<double> input;
  double value;
  while (std::cin >> value)  
    input.push_back(value);

  // for (int i = 0; i < input.size(); i++) {
  //   std::cout << input[i] << ", ";
  // }
  // std::cout << std::endl;
  // std::cout << "done" << std::endl;

  std::cout << "input size = " << input.size() << std::endl;

  // take in N
  int N = 16;
  std::cout << "please input a value for N (number of peaks): ";
  std::cin >> N;
  std::cout << "N is " << N << std::endl;

  CArray peaks(N * int(input.size()/WINDOW_LENGTH)); 
  int peak_index = 0; // this will go from 0 to 16. at each peaks[peak_index] will be a CArray of N sorted freq/amp vals

  /* for each sound clip of 2048 samples print information on the the N peaks.
     print the N frequency/amplitude pairs on a single line. separate the
     elements of each pair using the slash "/" character and separate each pair
     using comma "," character. For example:

     12.005/0.707,24.01/0.51,48.02/0.3 ...
     13.006/0.706,26.01/0.50,52.02/0.29 ...

     take N as an argument on the command line. -> done before while loop so that it is only done once
  */

  // grab every 2048 samples from input, overlapping 50%, until we get to the end of the file
  int input_index = 0; 
  int window_index = 0;
  double window[WINDOW_LENGTH]; //2048
  while ( input_index < input.size() ) {
    //std::cout << "window_index = " << window_index << std::endl;
    if (window_index >= WINDOW_LENGTH) { 
      window_index = 0; 
      input_index -= 0.5 * WINDOW_LENGTH; // 50% overlap

      //apply hanning window -> applies a sort of fading at the beginning and end of the window
      for (std::size_t i = 0; i < WINDOW_LENGTH; i++) { window[i] = hann(window[i]); }

      //zero pad the data
      CArray window_padded(PADDED_LENGTH); //8192
      for (std::size_t i = 0; i < PADDED_LENGTH; i++) {
        if (i < WINDOW_LENGTH) { window_padded[i] = window[i]; } else { window_padded[i] = 0; }
        //std::cout << window_padded[i] << std::endl;
      }
      
      //now, take an fft (forward)
      fft(window_padded);
      //std::cout << "fft" << std::endl;

      //let's get the peaks
      freq_amp pairs[PADDED_LENGTH];
      //this is hacky but working with valarray<Complex> crashes with sort
      for (auto& elem : window_padded) {
        //std::cout << elem << " ";
         // if the current index is needed:
        auto i = &elem - &window_padded[0];
        pairs[i] = std::make_pair(elem.real(), elem.imag());
        // any code including continue, break, return
      }
      
      std::sort(std::begin(pairs), std::end(pairs), std::greater<>()); // works now (yay)      

      // get the first N of these
      for (std::size_t i = 0; i < N; i++) {
        peaks[peak_index] = Complex(pairs[i].first, pairs[i].second);
        if (peak_index >= N) { peak_index += i + 1; } else { peak_index += i; } 
      }
    }
    window[window_index] = input[input_index];
    input_index++; window_index++;
  }
  print_data(peaks);
  return 0;
}