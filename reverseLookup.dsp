declare name "reverseLookup";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPLv3";

import("stdfaust.lib");



// TODO: special case duration = 0. run the algo from durationSlider=1 upwards;


// TODO: make tester: go trough all the parameter combinations and for every combination, increase the precision untill it hits the bell
// TODO: make 2 table.val for the sliders, and crosssfade them for the lookupVal, (or 4, for cub)
//
process =
  // reverseLookup(startFunInput, endFunInput, nrCompares, lookupFunc, lookupVal)
  // reverseLookupRaw(startFunInput, endFunInput, nrCompares, lookupFunc, lookupVal)
  // reverseLookupNdRaw(startFunInput, endFunInput, nrCompares, lookupFuncNd, shapeSlider, durationSlider, lookupVal)
  reverseLookupNd(startFunInput, endFunInput, nrCompares, lookupFuncNd, shapeSlider, durationSlider, lookupVal)
  : (_:hbargraph("pre-func", startFunInput, endFunInput))
    // : lookupFunc
  : lookupFuncNd(shapeSlider,durationSlider)
    <:(
    (
      (abs(_-lookupVal)
       // < (precision *ma.EPSILON)
       < (hslider("precision", 0, 0, 1, 0.001)*0.001)
      ):hbargraph("post = lookupFunc(pre)", 0, 1))
  , (_
     // :pow(0.1)
     :hbargraph("post-func", startFunOutput, endFunOutput))
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




reverseLookupNd(startFunInput, endFunInput, nrCompares, lookupFunc, shapeSlider, duration, lookupVal) =
  ba.tabulateNd(0,
                reverseLookupNdRaw(startFunInput, endFunInput, nrCompares, lookupFuncNd )
                , (Sy, Sdur, Sx, ry0, rDur0, startFunOutput, ry1, rDur1, endFunOutput, shapeSlider, duration, lookupVal)).val
with {
  // Sx = 1<<8;
  // Sx = 1<<12;
  // Sx = 1<<16;
  Sx = nrVals+1;
  // Sy = (nrShapes/shapeStep)+1;
  Sy = nrShapes+1;
  Sdur = nrDurations+1;
  ry0 = 0;
  rDur0 = 0;
  ry1 = nrShapes;
  rDur1 = nrDurations;
};

durationSlider = hslider("duration", 1, 0, nrDurations, 1);
dur2sec(x) = x/nrDurations:pow(2)*maxSeconds;
nrDurations = 32;
maxSeconds = 1;

reverseLookupNdRaw(startFunInput, endFunInput, nrCompares, lookupFuncNd, y, duration, lookupVal) =
  (startFunInput,endFunInput)
  : seq(i, nrCompares, compare)
  : +:_*.5 // average start and end
with {
  compare(start,end) =
    select2(bigger , start , middle)
  , select2(bigger , middle , end)
  with {
  bigger = lookupVal>lookupFuncNd(y,duration,middle);
  middle = (start+end)*.5;
};
};



endFunInput =
  1;
// 4096;
// 2;
lookupFunc(x) =
  (par(i, N,x/(i+1): warpedSineFormula(half)):>_/3)
;
N = 300;
lookupFuncNd(y,duration,x) =
  // warpedSineFormula(y,x);
  // shapeDifFormula(y,phase,duration,sr)
  // shapeDifFormula(y,x,duration:dur2sec,sr)
  rampCompare(y,x,duration:dur2sec,sr)
with {
  // phase = 0.5;
  sr = 48000;
  // sr = 1;
};

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
// lookupVal = hslider("lookupVal", startFunOutput, startFunOutput, endFunOutput, 0.001);
lookupVal = (hslider("lookupVal", 0, 0, nrVals, 1)/nrVals)
            *(endFunOutput-startFunOutput)
            +startFunOutput;
nrVals = 1000;
precision =
  // 100;
  hslider("precision", 5, 1, 1000, 1) * 1000000000;

sineShaper(x) = (sin((x*0.5 + 0.75)*2*ma.PI)+1)*0.5;
startFunInput = 0;
startFunOutput =
  0;
// lookupFuncNd(0,0.5,10);
// lookupFuncNd(0,nrDurations,startFunInput);
endFunOutput =
  // lookupFuncNd(0,0.5,9);
  1;
// 0.05:pow(0.5);
// faust2svg -double -sd reverseLookup.dsp && xdg-open reverseLookup-svg/process.svg
// turns out it's 1!
// par(i, nrShapes, par(j, nrDurations, par(k, nrPhases, lookupFuncNd(i,j:max(ma.EPSILON),k/nrPhases))))
// :ba.parallelMax(nrShapes*nrDurations*nrPhases) ;
// nrPhases = 48000*maxSeconds;
nrPhases = 16;
//lookupFunc(endFunInput);
// nrCompares = 22;
// nrCompares = 24;
nrCompares = 32;
// nrCompares = 64;
// nrCompares = 128;
// nrCompares = 256;
// nrCompares = 512;
// nrCompares = 1024;

nrShapes = 8;
half = nrShapes*.5;
shapeStep = 0.1;
shapeSlider = half+hslider("[2]release shape" , 0, 0-half, half, shapeStep);
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
  / nrShapes
  * range
  + start
  // : hbargraph("shapeBG", 0.3, 0.7)
with {
  range = 2* (.5-start);
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

// shapeDifFormula(shapeSlider,phase,0,sr) = 1;

shapeDifFormula(shapeSlider,phase,duration,sr) =
  warpedSineFormula(shapeSlider,phase+(1 / sr
                                       / duration
                                      )
                                // :min(1)
                   )
  - warpedSineFormula(shapeSlider,phase)
  // :min(1)
;

// warpedSineFormula(shapeSlider,phase)*(1+duration);
// warpedSineFormula(shapeSlider,phase)*(duration);
// phase*duration;

rampCompare(shapeSlider,phase,duration,sr) =
  shapeDifFormula(shapeSlider,phase,duration,sr)
  * (1/(1-warpedSineFormula(shapeSlider,phase)))
  // :min(1)
;
