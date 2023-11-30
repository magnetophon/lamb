declare name "reverseLookup";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPLv3";

import("stdfaust.lib");

process =
  reverseLookup(startFunInput, endFunInput, nrCompares, lookupFunc, lookupVal)
  : (_:hbargraph("pre-func", startFunInput, endFunInput))
  : lookupFunc
    <:(
    (
      (abs(_-lookupVal)
       < (precision *ma.EPSILON)
      ):hbargraph("post = lookupFunc(pre)", 0, 1))
  , (_:hbargraph("post-func", startFunOutput, endFunOutput))
  )
;
precision =
  // 100;
  hslider("precision", 4, 1, 100, 1);

// which x should I give to
// lookupFunc(x)
// so that the output is lookupVal
//
// lookupFunc should be rising only
reverseLookup(startFunInput, endFunInput, nrCompares, lookupFunc, lookupVal) =
  (startFunInput,endFunInput)
  : seq(i, nrCompares, compare)
  : +:_*.5 // average start and end
with {
  compare(start,end) =
    select2(bigger , start , middle)
  , select2(bigger , middle , end)
  with {
  bigger = lookupVal>lookupFunc(middle);
  middle = (start+end)*.5;
};
};

// endFunInput = 1;
// lookupFunc(x) = sineShaper(x^(1/div))^div;
endFunInput = div;
lookupFunc(x) = sineShaper(x/div);
// lookupFunc = _/div;
div =
  maxDiv;
// hslider("div", 1, 1, maxDiv, 1);
maxDiv = 4;
lookupVal = hslider("lookupVal", startFunOutput, startFunOutput, endFunOutput, 0.01);
sineShaper(x) = (sin((x*0.5 + 0.75)*2*ma.PI)+1)*0.5;
startFunInput = 0;
startFunOutput = lookupFunc(startFunInput);
endFunOutput = lookupFunc(endFunInput);
nrCompares = 32;
