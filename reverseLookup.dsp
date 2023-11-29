declare name "reverseLookup";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPLv3";

import("stdfaust.lib");

process =
  // os.lf_sawpos(0.5)
  lookupVal
  : reverseLookup(start, end, nrCompares, lookupFunc)
  : lookupFunc
    <:(
    ((abs(_-lookupVal)<(hslider("mult", 1, 1, 100, 1)*ma.EPSILON)):hbargraph("check", 0, 1))
  , hbargraph("val", start, end)
  , (lookupFunc(lookupVal):hbargraph("wrong", 0, 1))
  );

// which x should I give to
// lookupFunc(x)
// so that the output is lookupVal
//
// lookupFunc should be rising only
reverseLookup(start, end, nrCompares, lookupFunc, lookupVal) =
  (start,end,lookupVal)
  : seq(i, nrCompares, compare)
  : ((+:_*.5),!) // average start and end, throw away the rest
with {
  compare(start,end,lookupVal) =
    (
      select2(bigger , start , middle)
    , select2(bigger , middle , end)
    , lookupVal
    )
  with {
  bigger = lookupVal>lookupFunc(middle);
  middle = (start+end)*.5;
};
};

lookupVal = hslider("lookupVal", start, start, end, 0.01);
start = 0;
end = 1;
nrCompares = 32;
lookupFunc = sineShaper;
sineShaper(x) = (sin((x*0.5 + 0.75)*2*ma.PI)+1)*0.5;
