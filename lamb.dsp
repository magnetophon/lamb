declare name "lamb";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPLv3";

// TODO:
//
//  when calculating which part of the fullDif we already did, take into account the prev fullDiff and prevRamp
//
//  optimise: make attack and release mult factors, to be applied after the tables, and in reverse to the index of the tables
//
// att/rel discontinuities:
//
// alt solution:
//
// if (needToGoDownWithinReleaseTime)
//   then dontGoUp
//
// - check if need to go down within 2*attack
// - if yes:
//   - when dir rel == dir att(neg saw)
//   - use neg saw
//
//   check shape at neg ramp
//
//   check if better at att==rel
//
// optimize CPU:
//   make lookup dual SR
//   make lookup simpler:
//   - make 2 tables:
//     - one for (each?) shape (updated controll rate?)
//     - one for (prevSpeed, totalGR and shape) -> phase

// import("/home/bart/source/lamb/stdfaust.lib");
// import("/home/bart/source/faustlibraries/stdfaust.lib");
import("/nix/store/mljn5almsabrsw6mjb70g61688kc0rqj-faust-2.68.1/share/faust/stdfaust.lib");
// import("stdfaust.lib");

sinfun(x) = sin(pow(4*x/16,2));
sfx = hslider("sfx", 0, 0, 16, 1);
process =
  // ba.slidingMin(0,2)
  // , ba.slidingMin(1,2)
  // , ba.slidingMin(n,2)
  // , ba.slidingMin(0,4)
  // , ba.slidingMin(1,4)
  // , ba.slidingMin(2,4)
  // , ba.slidingMin(3,4)
  // , ba.slidingMin(n,4)
  // , ba.slidingMin(n,8)
  // ;
  // n = hslider("n", 0, 0, 8, 1);
  // sin(2*ma.PI*os.lf_saw(0.7));
  AR_tester;
// ba.tabulate(0, sinfun, 16, 0,16, sfx).lin;
// ,sinfun(sfx)
// ;
// co.FFcompressor_N_chan(strength,thresh,attack,release,knee,prePost,link,meter,2);
// co.RMS_FBFFcompressor_N_chan(strength,thresh,att,rel,RMStime,knee,prePost,link,FBFF,meter,2);
// co.RMS_FBcompressor_peak_limiter_N_chan(strength,thresh,threshLim,att,rel,RMStime,knee,link,meter,meterLim,2);

// par(i, 2, _*ba.db2linear(hslider("input gain", 0, -24, 24, 1):si.smoo)):
// lookahead_compressor_N_chan(strength,thresh,attack,release,knee,link,meter,2) ;

AR_tester =
  hgroup("",
         vgroup("[2]test", test)

         <:vgroup("[1]AR",
                  (AR(attack,release)
                   // : par(i, 3, _@(maxSampleRate-attackSamples))
                  )

                  // ,_@attackSamples
                  ,_@maxSampleRate
                   // ,ba.slidingMin(attackSamples,maxSampleRate)
                   // ,(((ba.slidingMin(attackSamples,maxSampleRate):smootherCascade(4, releaseOP, attackOP )),_@maxSampleRate):min)
                 ));
att = attack;
rel = release;
meterLim =meter;
threshLim = AB(threshLimP);
threshLimP = hslider("[03]thresh lim",0,-30,6,1);
RMStime = AB(RMStimeP);
RMStimeP = hslider("[03]RMS time",5,0,100,1)*0.001;
maxSampleRate = 192000;

AR(attack,release) =
  // (negative_ramp~_),

  loop~(_,_,_,_,_)
       :(_,_,_,!,!)
        // :(_,_@releaseSamples,_)
        // : (
        // par(i, 2, _@(maxSampleRate-attackSamples))
        // ,_)
        // : par(i, 3, _@(maxSampleRate-attackSamples))
with {
  loop(prevRamp,prevGain
       , prevHold
       ,prevAttacking,prevReleasing
       ,x) =
    ramp
  , gain
    // , attHold
    // , hold
  , limitMe
  , attacking
  , releasing
    // , negRamp
    // , warpNegRamp
    // , x
    // , (x:seq(i, 3, si.onePoleSwitching(releaseOP,attackOP)))
    // , (x==gain)
  with {

  attackHold = attack:ba.sAndH(1-prevAttacking);
  releaseHold = release:ba.sAndH(1-prevReleasing);
  holdHold = holdTime:ba.sAndH(1-prevReleasing);
  // attackHold =
  // attack:ba.sAndH(prevGain==(x@(maxSampleRate))) ;
  // releaseHold =
  // release:ba.sAndH(prevGain==(x@(maxSampleRate))) ;
  // attackSamples = (ba.sec2samp(attackHold)@_:ba.sAndH(1-prevAttacking))~(max(0,maxSampleRate-_):min(maxSampleRate));
  attackSamples = ba.sec2samp(attackHold);
  releaseSamples = ba.sec2samp(releaseHold);
  holdSamples =
    attackSamples +
    releaseSamples;
  // holdSamples = ba.sec2samp(holdHold);
  duration =
    // select3(attacking+releasing*2,1,attackHold,releaseHold);
    ((attackHold*attacking)+(releaseHold*releasing))
    // :max(ma.EPSILON)
  ;
  gain = prevGain+gainStep
         // :min(x@maxSampleRate)
  ;
  gainStep =
    select2(releasing
           , rawGainStep :max(dif)
           , rawGainStep :min(dif)
           ) with {
    rawGainStep =
      // select2(duration==0
      // ,
      shapeDif(shapeSlider,ramp,duration,ma.SR)*fullDif
      // , dif)
    ;
    fullDif =
      dif/(1-warpedSine(shapeSlider,ramp));
  };
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


  hold =
    // attHold
    fancyHold
  ;
  longHold =
    ba.slidingMin(holdSamples,maxSampleRate,x
                                            @max(0,(maxSampleRate-holdSamples)));

  switch = (prevGain>longHold) & (prevGain<=attHold);
  // switchStart =

  limitMe = prevGain<
            select2(checkbox("longHoldlim")
                   , attHold
                   ,longHold);
  fancyHold =
    longHold
    // <:select2(limitMe,_,max(prevGain))
    <:select2(limitMe,max(prevGain),_)
      // <:select2(checkbox("lim"),_,max(prevGain))
    :min(attHold)
  ;
  predictedGain= prevGain+((prevGain')-(prevGain''));
  // predictedGain= prevGain+((prevGain'')-(prevGain'));

  attHold = x
            @max(0,(maxSampleRate-attackSamples)):
            // @max(0,(maxSampleRate-holdSamples)):
            ba.slidingMin(attackSamples+1,maxSampleRate);

  // dif = attHold-prevGain;
  dif = hold-prevGain;
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
    :seq(i, 18, compare)
    : ((+:_*.5),!) // average start and end, throw away the rest
      // :max(start):min(end)
  with {
    start = 0;
    end = 1;
    rampStep = 1 / ma.SR / duration;
  };
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
    // SIZE = 1<<17; // for 2d lin
    SIZE = 1<<16;
    // table(i) = ba.tabulate(0, warpedSineFormula(shapeSliderVal(i)), SIZE, 0, 1, x).lin;
    // 16 compares: 4.5 -6%CPU
    // 21 compares: 7 % CPU
    // SIZE = 1<<9;
    // SIZE = 1<<3;
    // table(i) = ba.tabulate(0, warpedSineFormula(shapeSliderVal(i)), SIZE, 0, 1, x).cub;
  };


  warpedSineFormula(shapeSlider,x) =
    // sineShaper(warp(shape,knee,x)):pow(power)
    sineShaper(warp(shape,knee,x:max(-1):min(1))):pow(power)
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


lookahead_compressor_N_chan(strength,thresh,att,rel,knee,link,meter,N) =
  si.bus(N) <: si.bus(N*2):
  (
    (par(i,N,abs) : lookahead_compression_gain_N_chan_db(strength,thresh,att,rel,knee,link,N))
   ,si.bus(N)
  )
  : (ro.interleave(N,2) : par(i,N,(meter: ba.db2linear)*(_@maxSampleRate)))
;

lookahead_compression_gain_N_chan_db(strength,thresh,att,rel,knee,link,1) =
  lookahead_compression_gain_mono_db(strength,thresh,att,rel,knee);

lookahead_compression_gain_N_chan_db(strength,thresh,att,rel,knee,link,N) =
  si.bus(N)
  <: (si.bus(N),(ba.parallelMax(N) <: si.bus(N))) : ro.interleave(N,2) : par(i,N,(it.interpolate_linear(link)))
  : par(i,N,lookahead_compression_gain_mono_db(strength,thresh,att,rel,knee)) ;

lookahead_compression_gain_mono_db(strength,thresh,att,rel,knee) =
  ba.linear2db : gain_computer(strength,thresh,knee)
  : ba.db2linear:AR(attack,release)
                 // ,prevNegRamp,tmp
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
releaseP = hslider("[05]release",80,0,1000,1)*0.001;
releaseShape = AB(releaseShapeP);
releaseShapeP = half+hslider("[2]release shape" , 0, 0-half, half, 0.1);
holdTime = AB(holdTimeP);
holdTimeP = hslider("[06]holdTime",80,0,1000,1)*0.001;
offset = AB(offsetP);
// offsetP = hslider("offset", 0, -48, 0, 0.1):ba.db2linear;
offsetP = hslider("offset", 0, -1, 0, 0.01);
knee = AB(kneeP);
kneeP = hslider("[07]knee",2,0,30,1);
link = AB(linkP);
linkP = hslider("[08]link", 0, 0, 100, 1) *0.01;
FBFF = AB(FBFFP);
FBFFP = hslider ("[09]fb-ff",100,0,100,1) *0.01;
power = AB(powerP);
powerP = hslider("[10]power", 2, ma.EPSILON, 10, 0.1);
PMItime = AB(PMItimeP);
PMItimeP = hslider("[11]PMI time",20,0,1000,1)*0.001;
dw = AB(dwP);
dwP = hslider ("[12]dry/wet",100,0,100,1) * 0.01:si.smoo;

attackOP = AB(attackOpP);
attackOpP = hslider("[13]attack OP",6,0,1000,1)*0.001;
releaseOP = AB(releaseOpP);
releaseOpP = hslider("[14]release OP",80,0,1000,1)*0.001;
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
