#!/bin/bash

function usage {
  echo -e "Usage:"
  echo -e "\t-i s: path to the raw PCM file"
  echo -e "\t-p s: path to output spectrogram image file"
  echo -e "\t-d s: path to ttf data text file"
  echo -e "\t-z d: number of Herz per bin in frequency domain, default: 5"
  echo -e "\t-b d: number of bins to reserve in frequency domain, default: 200"
  echo -e "\t-f d: feature flag"
  echo -e "\t-h: print this help info"
  exit 1
}

audio_file=
specgrm_file=
fft_data_file=
herz_per_bin=
bins_to_reserve=
round_to_avg=
feature_flag=-1

while getopts "i:p:d:z:b:r:f:h" opt "$@"; do
  case $opt in
    i) 
      audio_file=$OPTARG; 
      ;;
    p) 
      specgrm_file=$OPTARG; 
      ;;
    d) 
      fft_data_file=$OPTARG; 
      ;;
    z) 
      herz_per_bin=$OPTARG; 
      ;;
    b) 
      bins_to_reserve=$OPTARG; 
      ;;
    r) 
      round_to_avg=$OPTARG; 
      ;;
    f) 
      feature_flag=$OPTARG; 
      ;;
    \?) 
      usage 
      ;;
    :) 
      usage
      ;;
  esac
done

if [[ ! -f "$audio_file" ]]; then
  echo "$audio_file does not exist"
  usage
fi

if [[ "$specgrm_file" = "" ]]; then
  echo "ERROR: must specify path to output spectrogram image file"
  usage
fi

if [[ "$fft_data_file" = "" ]]; then
  echo "ERROR: must specify path to output fft data file"
  usage
fi

if [[ "$feature_flag" -eq -1 ]]; then
  echo "ERROR: must specify a feature flag, 0 or 1"
  usage
fi

tmp_raw_file=$(echo $audio_file | awk -F. '{for(i=1;i<NF;++i){printf $i"."}}')raw

# read sample_rate and number of channels
read sample_rate is_stereo <<< $(ffmpeg -i "$audio_file" 2>&1 | perl -ne 'if(/(\d+) Hz/){print $1; if(/2 channels|stereo/){print " 1\n"}else{print " 0"}}')

if [[ ! -f "$tmp_raw_file" ]]; then
  # convert the audio file into raw PCM
  echo "converting $audio_file into raw PCM..."
  ffmpeg -i "$audio_file" -f s16le -acodec pcm_s16le "$tmp_raw_file" 2> /dev/null
fi


echo "Processing..."

./aud2spctrogram -i "$tmp_raw_file" -p "$specgrm_file" -d "$fft_data_file" -z $herz_per_bin -b $bins_to_reserve -s $sample_rate -c $is_stereo -r $round_to_avg -f $feature_flag

echo "     Input: $audio_file"
echo "Spctrogram: $specgrm_file"
echo "  FFT data: $fft_data_file"
echo "SampleRate: $sample_rate"
echo "Is Stereo?: $is_stereo"
echo
