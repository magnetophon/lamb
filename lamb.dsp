declare name "lamb";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPLv3";

import("/home/bart/source/lamb/stdfaust.lib");

sinfun(x) = sin(pow(4*x/16,2));
sfx = hslider("sfx", 0, 0, 16, 1);
process =
  AR_tester;
// par(i, 2, _*ba.db2linear(hslider("input gain", 0, -24, 24, 1):si.smoo)):
// lookahead_compressor_N_chan(strength,thresh,attack,release,knee,link,FBFF,meter,2);

AR_tester =
  hgroup("",
         vgroup("[2]test", test)
         <:vgroup("[1]AR",
                  (AR(attack,release) :(!,_,_,_))
                  ,_@maxHoldTime
                   // ,ba.slidingMin(attackSamples,maxSampleRate)
                   // ,(((ba.slidingMin(attackSamples,maxSampleRate):smootherCascade(4, releaseOP, attackOP )),_@attackSamples):min)
                 ));
att = attack;
rel = release;
attackSamples = ba.sec2samp(attack);
releaseSamples = ba.sec2samp(release);
meterLim =meter;
threshLim = AB(threshLimP);
threshLimP = hslider("[03]thresh lim",0,-30,6,1);
RMStime = AB(RMStimeP);
RMStimeP = hslider("[03]RMS time",5,0,100,1)*0.001;
maxSampleRate = 192000;
maxHoldTime = 2 * maxSampleRate;

AR(attack,release) = loop~(_,_,_)
                          // :(!,_)
with {
  loop(prevRamp, prevTurnAroundRamp, prevGain, x) =
    ramp
  , turnAroundRamp
  , turnAroundGain
    // , gain
    // , turnAroundHold
  , (turnAroundGainStep:(!,_))
    // , attackHold
    // , downSoon
    // , select2(speedMatch,gain,turnAroundGain)
    // , turnAroundGain
    // , x
    // , (x:seq(i, 3, si.onePoleSwitching(releaseOP,attackOP)))
    // , (x==gain)
  with {
  duration =
    // select3(attacking+releasing*2,1,attack,release);
    (attack*attacking)+(release*releasing);
  gain = prevGain+gainStep ;
  gainStep =
    select2(releasing
           , rawGainStep :max(dif)
           , rawGainStep :min(dif)
           ) with {
    rawGainStep =
      shapeDif(shapeSlider,ramp,duration,ma.SR)*fullDif;
  };
  fullDif =dif/(1-warpedSine(shapeSlider,ramp));
  shapeDifFormula(shapeSlider,phase,len) =
    warpedSineFormula(shapeSlider,phase+len)
    - warpedSineFormula(shapeSlider,phase);

  shapeDif(shape,phase,duration,sr) =
    // ba.tabulateNd(1,shapeDifFormula,(nrShapes,1<<17,1<<7,0,0,1/maxSampleRate/1,nrShapes,1,1/24000/(1/maxSampleRate),shapeSlider,phase,(1 / ma.SR / duration))).lin;
    // ba.tabulateNd(0,shapeDifFormula,(nrShapes,1<<17,1<<7,0,0,1/maxSampleRate/1,nrShapes,1,1/24000/(1/maxSampleRate),shapeSlider,phase,(1 / ma.SR / duration))).lin;
    // ba.tabulateNd(0,shapeDifFormula,(nrShapes,1<<17,1<<7,0,0,1/maxSampleRate/1,nrShapes,1,1/24000/(1/maxSampleRate),shapeSlider,phase,(1 / ma.SR / duration))).lin;
    // ba.tabulateNd(1,shapeDifFormula,(3,1<<16,1<<6,0,0,1/48000/1,nrShapes,1,1/24000/(1/48000),shapeSlider,phase,(1 / ma.SR / duration))).lin;
    warpedSine(shapeSlider,phase+(1 / sr / duration))
    - warpedSine(shapeSlider,phase);
  // warpedSineFormula(shapeSlider,phase+(1 / sr / duration))
  // - warpedSineFormula(shapeSlider,phase);
  attackHold =
    ba.slidingMin(     attackSamples,      maxSampleRate,x)@(maxHoldTime-(attackSamples));
  releaseHold =
    ba.slidingMin(     releaseSamples,      maxSampleRate,x)@(maxHoldTime-(releaseSamples));
  turnAroundHold =
    releaseHold;
  // max(attackHold,releaseHold);
  dif = attackHold-prevGain;
  releasing =
    dif>0;
  attacking =
    dif<0;

  compare(start,end,compSlope) =
    (
      select2(bigger , start , middle)
    , select2(bigger , middle , end)
    , compSlope
    )
  with {
    bigger = compSlope>slope(middle);
    slope(x) =
      shapeDif(shapeSlider,x,duration,ma.SR)
      *(1/(1-warpedSine(shapeSlider,x)));
    middle = (start+end)*.5;
  };
  // test with shape minimal, so 0.3 and duration = (3/16)^2
  // (lower shapes give jumps in the phase anyway)
  // at 48k, 13 seems to little, 14 works
  //
  // test with shape -0.4, and duration = (10/16)^2
  // at 48k, 14 seems to little, 15 works
  //
  // test with shape 3.2, and duration = (1/16)^2
  // at 48k, 16 seems to little, 24 works
  //
  // 15 takes about as much CPU as 16, so better be safe than sorry for now
  //
  // at 406.5 ms, we get a too slow ramp with 18 compares
  // 20 is ok, 22 is closer, 21 is "good enough"TM and cheaper
  //
  // with the above settings, too low nr of compares gives a stuck or too slow ramp
  ramp =
    (start,end)
  , shapeDif(shapeSlider,prevRamp+rampStep,duration',ma.SR)
    * ((dif'/dif)/(1-warpedSine(shapeSlider',prevRamp)))
    :seq(i, 22, compare)
    : ((+:_*.5),!) // average start and end, throw away the rest
    :max(0):min(1)
  with {
    start = 0;
    end = 1;
    rampStep = 1 / ma.SR / duration;
  };
  turnAroundRamp = startTurnAroundRamp *
                   ( prevTurnAroundRamp + turnAroundRampStep);
  // with {
  startTurnAroundRamp = goingUp & downSoon;
  goingUp = gain > prevGain;
  downSoon =
    (turnAroundHold < attackHold)
    // & (turnAroundHold < prevGain)
    // & (turnAroundHold<gain)
  ;

  turnAroundRampStep = 1 / ma.SR / attack;
  // };

  // speedMatch = turnAroundSpeed >= (gain - prevGain);
  speedMatch =
    (turnAroundGain >= prevGain)
    & startTurnAroundRamp;
  turnAroundGain =
    (prevGain + (turnAroundGainStep:(_,!)))
    // : smoothAtZero
  ;
  smoothAtZero = si.onePoleSwitching(swRel,0);
  swRel = turnAroundRamp >= 1;
  turnAroundGainStep = loop~_:(!,_,_)
  with {
    loop(prevStartValue) =
      startValue
    , (gainStep * mult)
    , mult
    with {
    startValue = select2(switchOn
                        , prevStartValue
                        , turnAroundRamp);
    switchOn = (turnAroundHold<prevGain) & ((turnAroundHold<prevGain)'==0);
    keepOn = (turnAroundRamp < 1) & (turnAroundRamp > 0);
    switch = select2(switchOn
                    , _*keepOn
                    , 1)~_;
    rawMult = 1-(switch*(
                    (turnAroundRamp-startValue)/(1-startValue)
                  ));
    mult = loop~_
    with {
      loop(prevMult) = select2(sel
                              , rawMult
                                // , (1-1'):max(7*turnAroundRampStep)
                              , (attacking):max(7*turnAroundRampStep)
                              )
      with {
      sel = (prevMult<(8*turnAroundRampStep)) & (prevGain<attackHold);
    };

    };

  };
  };

  // warpedSine(shapeSlider,phase+(1 / sr / duration))
  // - warpedSine(shapeSlider,phase);
  // (warpedSine(attackShape,(turnAroundRamp)+(1 / ma.SR / duration))
  // - (warpedSine(attackShape,(turnAroundRamp))))
  // * fullTurnAroundDif;
  // * fullDif;
  fullTurnAroundDif =turnAroundDif/(1-warpedSine(attackShape,turnAroundRamp));
  // fullTurnAroundDif =turnAroundDif/(1-warpedSine(attackShape,ramp));

  turnAroundDif = turnAroundHold-prevGain;


  // ******************************************** the curves: ******************************
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
  warpedSine(shapeSlider,x) =
    // at low number of compares the raw formula is faster than the tabulated version
    // 16 compares: 5 to 6 % CPU
    // when doing contant phase recalculations, we need higher precision in the newramp function
    // cause we get wrong ramp durations (to steep or not steep enough) otherwise
    // 21 compares seems to work well enough in all cases so far
    // at the higher number of compares (21) we get 11-12% CPU for the raw formaula
    // warpedSineFormula(shapeSlider,x)
    // the tables do much better
    // Size can be 1<<3;
    // ba.tabulateNd(1, warpedSineFormula,(nrShapes, 1<<3,0, 0,nrShapes, 1, shapeSlider,x)).cub
    ba.tabulateNd(0, warpedSineFormula,(nrShapes, SIZE,0, 0,nrShapes, 1, shapeSlider,x)).lin
    // par(i, nrShapes+1, table(i) * xfadeSelector(shapeSlider,i)):>_
    // this one is only slightly cheaper, but less user freindly
    // par(i, nrShapes+1, table(i) * ((shapeSlider)==i)):>_
  with {
    // with 16 compares: 4.5 to 5.5 % CPU
    // with 21 compares: 12 - 17 % CPU!
    // 23 goes wrong with rel=3, shape=minimum, ramp not steep enough
    // SIZE =24 hangs the ramp with rel=406.5 ms, any shape
    // SIZE>24 doesn't compile, even with -quad
    // SIZE = 1<<24;
    // table(i) = ba.tabulate(0, warpedSineFormula(shapeSliderVal(i)), SIZE, 0, 1, x).val;

    // 16 compares: 3 to 4 % CPU
    // 21 compares: 4.5 to 5.5 % CPU
    // test with
    // patho case: rel 1 shape -3.4
    // patho case: rel 1 shape -1.8
    SIZE = 1<<17; // for 2d lin
    // SIZE = 1<<16;
    // table(i) = ba.tabulate(0, warpedSineFormula(shapeSliderVal(i)), SIZE, 0, 1, x).lin;
    // 16 compares: 4.5 -6%CPU
    // 21 compares: 7 % CPU
    // SIZE = 1<<9;
    // SIZE = 1<<3;
    // table(i) = ba.tabulate(0, warpedSineFormula(shapeSliderVal(i)), SIZE, 0, 1, x).cub;
  };


  warpedSineFormula(shapeSlider,x) =
    // sineShaper(warp(shape,knee,x)):pow(power)
    sineShaper(warp(shape,knee,x:max(0):min(1))):pow(power)
  with {
    power = (4*shape/3)+(1/3);
    knee = min(2*shape,2-(2*shape));
    shape = shapeSliderVal(shapeSlider);
  };
  shapeSlider =
    // select2(releasing, 1-slider)
    select2(releasing
           , attackShape
           , releaseShape);


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

  ishape =
    shapeSliderVal(shapeSlider);
};
};


lookahead_compressor_N_chan(strength,thresh,att,rel,knee,link,FBFF,meter,N) =
  si.bus(N) <: si.bus(N*2):
  (
    ((ro.interleave(N,2) : par(i,N*2,abs) :par(i,N,it.interpolate_linear(FBFF))
      : lookahead_compression_gain_N_chan_db(strength*(1+((FBFF*-1)+1)),thresh,att,rel,knee,link,N)),si.bus(N))
    : (ro.interleave(N,2) : par(i,N,(meter: ba.db2linear)*(_@attackSamples)))
  )
  ~ si.bus(N);

lookahead_compression_gain_N_chan_db(strength,thresh,att,rel,knee,link,1) =
  lookahead_compression_gain_mono_db(strength,thresh,att,rel,knee);

lookahead_compression_gain_N_chan_db(strength,thresh,att,rel,knee,link,N) =
  par(i,N,lookahead_compression_gain_mono_db(strength,thresh,att,rel,knee))
  <: (si.bus(N),(ba.parallelMin(N) <: si.bus(N))) : ro.interleave(N,2) : par(i,N,(it.interpolate_linear(link)));

lookahead_compression_gain_mono_db(strength,thresh,att,rel,knee) =
  ba.linear2db : gain_computer(strength,thresh,knee)
  : ba.slidingMin(attackSamples,maxSampleRate)
  : ba.db2linear:AR(attack,release)
  :(!,_) // for testing
  :ba.linear2db
with {
  gain_computer(strength,thresh,knee,level) =
    select3((level>(thresh-(knee/2)))+(level>(thresh+(knee/2))),
            0,
            ((level-thresh+(knee/2)) : pow(2)/(2*max(ma.EPSILON,knee))),
            (level-thresh))
    : max(0)*-strength;
};

AB(p) = ab:hgroup("[1]A/B",sel(aG(p),bG(p)));
sel(a,b,x) = select2(x,a,b);
aG(x) = vgroup("[0]a", x);
bG(x) = vgroup("[1]b", x);

B = si.bus(2);
ab = checkbox("[0]a/b");
bypass = AB(bypassP);
bypassP = checkbox("[00]bypass");
prePost = AB(prePostP);
prePostP = checkbox("[01]prePost");
strength = AB(strengthP);
strengthP = hslider("[02]strength", 100, 0, 100, 1) * 0.01;
thresh = AB(threshP);
threshP = hslider("[03]thresh",0,-30,6,1);
attack = AB(attackP);
attackP = hslider("[04]attack",6,0,1000,1)*0.001;
attackShape = AB(attackShapeP);
attackShapeP = half+hslider("[2]attack shape" , 0, 0-half, half, 0.1);
// release = AB(releaseP);
release = AB(releaseP);
releaseP = hslider("[05]release",80,1,1000,1)*0.001;
releaseShape = AB(releaseShapeP);
releaseShapeP = half+hslider("[2]release shape" , 0, 0-half, half, 0.1);
knee = AB(kneeP);
kneeP = hslider("[06]knee",2,0,30,1);
link = AB(linkP);
linkP = hslider("[07]link", 0, 0, 100, 1) *0.01;
FBFF = AB(FBFFP);
FBFFP = hslider ("[08]fb-ff",100,0,100,1) *0.01;
power = AB(powerP);
powerP = hslider("[09]power", 2, ma.EPSILON, 10, 0.1);
PMItime = AB(PMItimeP);
PMItimeP = hslider("[10]PMI time",20,0,1000,1)*0.001;
dw = AB(dwP);
dwP = hslider ("[11]dry/wet",100,0,100,1) * 0.01:si.smoo;

attackOP = AB(attackOpP);
attackOpP = hslider("[12]attack OP",6,0,1000,1)*0.001;
releaseOP = AB(releaseOpP);
releaseOpP = hslider("[12]release OP",80,0,1000,1)*0.001;
nrShapes = 9;
half = (nrShapes-1)*.5;

ARtest = toggle(soft,loud) with {
  toggle(a,b) = select2(block,b,a);
  block = os.lf_sawpos(0.5)>0.5;
  soft = sine*0.1;
  loud = sine;
  sine = os.osc(5000);
};


meter =
  _<: attach(_, (max(-24):min(0):hbargraph(
                   "v:[10]meters/[unit:dB]", -24, 0)
                ));

///////////////////////////////////////////////////////////////////////////////
//                                    test                                   //
///////////////////////////////////////////////////////////////////////////////

test = (select3(hslider("test", 2, 0, 2, 1)
               , test0
               , test1
               , test2
               )

       , no.lfnoise(hslider("rate", 100, 0.1, 20000, 0.1))
       )
       :it.interpolate_linear(hslider("Xfade", 0, 0, 1, 0.01))
;

test0 = select2(os.lf_sawpos(0.5)>0.5, -1,1);
test1 = select3(
          ((os.lf_sawpos(1)>hslider("POS1", 0.25, 0, 1 , 0.01))+(os.lf_sawpos(1)>hslider("POS2", 0.5, 0, 1 , 0.01))),
          1, -1, 0);
test2 =
  (loop~_)
with {
  loop(prev,x) = no.lfnoise0(abs(prev*69)%9:pow(2)+1);
};
// TODO:
// will algo
// variable rms size
// rms as attack
// shaper around attack/release
// https://www.desmos.com/calculator/xeystvebfz
// https://www.desmos.com/calculator/ir2xakmtav
// sine shaper: https://www.desmos.com/calculator/a9esb5hpzu
// mult by shaper to get gain, div by shaper to calc dif
// piecewise:
// https://www.desmos.com/calculator/kinlygqpcf
// knee:
// https://www.desmos.com/calculator/b3nqkydhye
// https://www.desmos.com/calculator/v6natnftsw
// https://www.desmos.com/calculator/zfksmqczif
// https://www.desmos.com/calculator/sbsdegqezh
// https://www.desmos.com/calculator/pobazzinqv
// https://www.desmos.com/calculator/znpkazx5zl
// https://www.desmos.com/calculator/jgi4uhuifs
// https://www.desmos.com/calculator/0m5ks4hraa
//
// complete:
// https://www.desmos.com/calculator/otuch9nfsc
// s(s(x)) with auto knee
// https://www.desmos.com/calculator/hlushh63kc
// both, with colors:
// https://www.desmos.com/calculator/u5jhxh5tk1
// compare to sine:
// https://www.desmos.com/calculator/zybaewqrai
// extra scaling:
// https://www.desmos.com/calculator/qk8gymjfsl
// scaling plus area under curve:
// https://www.desmos.com/calculator/vcl8vrty1yi
// no scaling area under curve
// https://www.desmos.com/calculator/wjtyrllnhd
// https://www.desmos.com/calculator/apeaxg6yxm
// add negative phase:
// https://www.desmos.com/calculator/nbf4dbuuj5
// derivative/speedDif:
// https://www.desmos.com/calculator/j7wk2oedsi

// TODO: stop ramp if we are not there yet on the steepest point.
// steepest => derivative of the derivative approaches 0.
// not there yet =>
// fullDif
//
// TODO:
// for when the max slope is not big enough:
// make a table of shape in to ramp at maxSlope out
// find (binary search) the shape that gives the wanted slope at the maxslope of that shape
// set ramp to that maxslope, shape offset untill done
// *****  OR  ******
// find the stepsize at which the slope matches at ramp=0.5

// TODO: continiously variable shape: see if we are going to make it and if not adapt shape, each sample
// TODO: when changerate too big, set shape to 0.5 and try again
// TODO: fix too slow speed at the beginning of short duration ramps when ramp is near ramp', but not equal: make a normal step.
// TODO: when ramp is zero, and gain<x : fade to x
// TODO: if you make the number of shapes the user can select small, say 16, you can use 16 lookup tables for the phase corrector
// TODO: use negative ramps when needed?
// for example when speed doesn't match up otherwise,
//   (this one needs an shape that always goes up)
// or when we're doing the release, and for the attack we need to change direction
//   (needs a shape that folows the sine, so at negative phases we have negative speed)
// TODO: for the shape difs at the outer edges, where it goes out of scope, use the values at the edges
// TODO: turnaround:
// when attacking and the going into release, keep attack ramp going untill speed is 0, then switch to release ramp
// TODO: makeup gain: implement as an offset to the wanted GR, before the smoothing, that way any automation is smoothet by us.
// same with strength
// TODO: link: before smoother
// TODO: binary search as a function lin the libraries
// TODO: auto makup gain by area under curve
// TODO: make sure we use ints where we can
// TODO: fix the out of bound reads from tabulateNd.cub:
// make the parameter write ranges a bit bigger than the read ranges
// TODO: much smaller number of compares, and if the error is small, just use the simple ramp
N=4;
T = ma.T;
PI = ma.PI;
TWOPI = 2.0 * PI;
TWOPIT = TWOPI * T;
/* Cascaded one-pole smoothers with attack and release times. */
smoother(N, att, rel, x) = loop ~ _
with {
  loop(fb) = coeff * fb + (1.0 - coeff) * x
  with {
  cutoffCorrection = 1.0 / sqrt(pow(2.0, 1.0 / N) - 1.0);
  coeff = ba.if(x > fb, attCoeff, relCoeff);
  TWOPITC = TWOPIT * cutoffCorrection;
  attCoeff = exp(-TWOPITC / att);
  relCoeff = exp(-TWOPITC / rel);
};
};
smootherCascade(N, att, rel, x) = x : seq(i, N, smoother(N, att, rel));

mysel(x)=         (checkbox("AR")*
                   (si.onePoleSwitching(hslider("rel simple", 8, 0, 1000, 0.1)*0.001
                                       ,hslider("att simple", 8, 0, 1000, 0.1)*0.001,x))
                  , (1-checkbox("AR"))*AR(x):(!,_,!)):>_,x;
