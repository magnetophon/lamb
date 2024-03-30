
declare name "lamb-rs";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPLv3";


// This is a slight adaption from lamb.dsp for usage in
// https://github.com/magnetophon/lamb-rs


process =
  component("lamb.dsp")
  [
    enableLatencyMeter = 1;
    // set this to 2 for downsampled selectable outputs,
    // one for level and one for gain reduction
    selectOutputs = 0;
  ];
