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

int main(int argc, char *argv[]) {
  // take the data in
  //
  vector<double> data;
  double value;
  int n = 0;
  while (std::cin >> value) {
    data[n] = value;
    n++;
  }

  // put your code here!
  //
  // for each sound clip of 2048 samples print information on the the N peaks.
  // print the N frequency/amplitude pairs on a single line. separate the
  // elements of each pair using the slash "/" character and separate each pair
  // using comma "," character. For example:
  //
  // 12.005/0.707,24.01/0.51,48.02/0.3 ...
  // 13.006/0.706,26.01/0.50,52.02/0.29 ...
  //
  // take N as an argument on the command line.
  //

  return 0;
}
