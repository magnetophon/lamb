
declare name "lamb-rs";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPLv3";


// This is a slight adaption from lamb.dsp for usage in
// https://github.com/magnetophon/lamb-rs


process =
  component("lamb.dsp")[enableLatencyMeter = 1;];
