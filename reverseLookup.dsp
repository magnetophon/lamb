declare name "reverseLookup";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPLv3";

import("stdfaust.lib");

process =
  reverseLookup(startFunInput, endFunInput, nrCompares, lookupFunc, lookupVal)
  // reverseLookupRaw(startFunInput, endFunInput, nrCompares, lookupFunc, lookupVal)
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
// which x should I give to
// lookupFunc(x)
// so that the output is lookupVal
//
// lookupFunc should be rising only
reverseLookup(startFunInput, endFunInput, nrCompares, lookupFunc, lookupVal) =
  // reverseLookupRaw(startFunInput, endFunInput, nrCompares, lookupFunc, lookupVal)
  ba.tabulate(0,
              reverseLookupRaw(startFunInput, endFunInput, nrCompares, lookupFunc)
              , S, startFunOutput, endFunOutput, lookupVal).lin
with {
  S = 1<<8;
  // S = 1<<12;
  // S = 1<<23;
};


reverseLookupRaw(startFunInput, endFunInput, nrCompares, lookupFunc, lookupVal) =
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

endFunInput = 1;
// lookupFunc(x) = sineShaper(x^(1/div))^div;
lookupFunc(x) =
  // (0.5*warpedSineFormula(shapeSlider,x))
  // +(0.5*
  // (x:seq(i, 3, warpedSineFormula(shapeSlider)))
  (par(i, N,x/(i+1): warpedSineFormula(shapeSlider)):>_/3)
  // )
;
N = 300;
// endFunInput = div;
// lookupFunc = _/div;
div =
  maxDiv;
// hslider("div", 1, 1, maxDiv, 1);
maxDiv = 4;
// patho case:
// tabulated with S = 1<<8; lookupVal=0.01 precision=228      0.5% - 0.8% CPU 21MiB
// tabulated with S = 1<<12; lookupVal=0.74 precision=145      0.7% - 0.8% CPU 1.7GB
// raw: lookupVal=0.68 precision=147   3.5% - 4% CPU 22MiB
lookupVal = hslider("lookupVal", startFunOutput, startFunOutput, endFunOutput, 0.01);
precision =
  // 100;
  hslider("precision", 8, 1, 1000, 1);

sineShaper(x) = (sin((x*0.5 + 0.75)*2*ma.PI)+1)*0.5;
startFunInput = 0;
startFunOutput = lookupFunc(startFunInput);
endFunOutput = lookupFunc(endFunInput);
// nrCompares = 24;
// nrCompares = 32;
nrCompares = 64;
// nrCompares = 128;
// nrCompares = 256;

nrShapes = 9;
half = (nrShapes-1)*.5;
shapeSlider = half;
warpedSineFormula(shapeSlider,x) =
  // sineShaper(warp(shape,knee,x)):pow(power)
  sineShaper(warp(shape
                 ,knee,x:max(-1):min(1))):pow(power)
with {
  power = (4*shape/3)+(1/3);
  knee = min(2*shape,2-(2*shape));
  shape = shapeSliderVal(shapeSlider);
};
// shapeSlider =
// select2(releasing
// , attackShape
// , releaseShape);


shapeSliderVal(shapeSlider) =
  shapeSlider
  / (nrShapes-1)
  * range
  + start
  // : hbargraph("shapeBG", 0.3, 0.7)
with {
  range = 2* (.5-start);
  // lower shapes then 0.3 give jumps in the phase at low durations (d < (3/16:pow(2)))
  // also they give stuck ramps at nr of compares < 14
  // shapeSliderVal(shapeSlider) = hslider("shape", 0.5, 0.30, 0.70, 0.01);
  start = 0.3;
};
// };
// };

kneeCurve(shape,knee,x) =
  select3( (x>shape-(knee*.5)) + (x>shape+(knee*.5))
         , 0
         , (x-shape + (knee*.5)):pow(2)/(knee*2)
         , x-shape);
warp(shape,knee,x) =
  (x-factor*kneeCurve(shape,knee,x))/(2*shape) with {
  factor = (1/shape-2)/(1/shape-1);
};
sineShaper(x) = (sin((x*0.5 + 0.75)*2*ma.PI)+1)*0.5;
