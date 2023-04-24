declare name "lamb";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPLv3";

import("stdfaust.lib");

// make a release based on:
// place

// crossfade between multiple inputs (just the gain)
process =
  // xfadeSelector(0.33,5)
  test:AR
       // test
       // PMI_FBFFcompressor_N_chan(strength,thresh,att,rel,knee,prePost,link,FBFF,meter,N);
       // (ARtest:PMI_compression_gain_mono_db(strength,thresh,att,rel,knee,prePost):ba.db2linear)
       // , os.lf_sawpos(1)>0.5
;


test0 = select3(
          ((os.lf_sawpos(1)>0.3)+(os.lf_sawpos(1)>0.5)),
          0, 1, -0.1);
test1 = select2(os.lf_sawpos(1)>0.5, 0.1,0.9);
test =
  ((loop~_)
  , no.lfnoise(hslider("rate", 100, 0.1, 1000, 0.1))
  )
  :it.interpolate_linear(hslider("Xfade", 0, 0, 1, 0.01))
with {
  loop(prev,x) = no.lfnoise0(abs(prev*69)%9:pow(2)+1);
};


AR = loop~(_,_)
          // :(!,si.bus(4))
with {
  loop(prevRamp,prevGain,x) =
    ramp
  , trueGain
    // rampShapeFix
    // , gainShapeFix
  , x
    // , warpedSine(shape,rawRamp)
    // , (gain==x)
    // , (gain==gain')
    // , prevFullDif
    // , coockedDif
    // , running
    // , releasing
    // , (abs(rawGainStep)/2)
    // , switchedShape
    // , newShape'
    // , shapeIntervention
    // , shapeInterventionHold
    // , newRampShapeX
    // ,(compareShape(0,1,0.5):>_)
    // , abs(changeRate)
    // ,(changeRate* (warpedSine(shape,rawRamp+rampStep)/(warpedSine(shape,rawRamp):max(smallest))))
    // ,(maxCR/ma.SR)
  , intervention
    // , (maxDerTable(shape) :hbargraph("MD", 0, 1))
    // , maxDerivative(ba.time/(1<<16))
  with {
  rawRamp = (prevRamp+rampStep)*running:min(1):max(0);
  ramp = select2(intervention,rawRamp,newRamp):max(0):min(1)*runningNotDirChange;
  rampStep = 1 / ma.SR / duration;
  duration = select3(attacking+releasing*2,1,attack,release);
  // attack = hslider("attack[scale:log]", 0.08, 0, 1, 0.001);
  attack =
    hslider("attack", 0.5, 0, nrRleases, 1 )/nrRleases:pow(2);
  // release = hslider("release[scale:log]", 0.5, 0.001, 1, 0.001):max(1 / ma.SR);
  release =
    hslider("release", 0.5, 0, nrRleases, 1 )/nrRleases:pow(2);
  nrRleases = 20;
  smallest = 1/192000;
  gain = prevGain+gainStep
         :max(-1):min(1);
  trueGain =
    select2( ((1-running)
              | (
                ((attack==0) & (rawDif<0))
                | ((release==0) & (rawDif>0))))
             & (1-dirChange)
           , gain
           , x
           );
  rawGainStep =
    shapeDif(shape,rawRamp,rampStep)*fullDif;
  shapeDif(shape,phase,step) =
    warpedSine(shape,phase+step)
    - warpedSine(shape,phase);

  gainStep = select2(releasing
                    , rawGainStep
                      :min(smallest)
                      :max(rawDif)
                    , rawGainStep
                      :max(smallest)
                      :min(rawDif)
                    )
             * running ;

  gainShapeFix = prevGain+gainStepShapeFix:max(-1):min(1);
  rawGainStepShapeFix =
    shapeDif(switchedShape,rawRamp,rampStep)*fullDifShapeFix;
  gainStepShapeFix = select2(releasing
                            , rawGainStepShapeFix:min(0-smallest):max(rawDif)
                            , rawGainStepShapeFix:max(smallest):min(rawDif)
                            )
                     * running;
  rawDif = x-prevGain;
  fullDif =rawDif/(1-warpedSine(shape,rawRamp));
  fullDifShapeFix =rawDif/(1-warpedSine(switchedShape,rawRamp));
  runningNotDirChange = running* (1-dirChange);
  running = (attacking | releasing) ;
  dirChange = (attacking != attacking')| (releasing != releasing');
  prevFullDif =rawDif/(1-warpedSine(shapeSliderVal(shapeSlider),prevRamp));
  coockedDif = (prevFullDif/ ((abs(prevGain-prevGain'):max(ma.EPSILON)*ma.SR)) );
  closeEnough = (abs(coockedDif)<=(5000*ma.EPSILON));
  // prevFullDif =rawDif/shape;
  releasing =
    prevGain<x;
  // (coockedDif>closeEnough);
  attacking =
    prevGain>x;
  // (coockedDif<(0-closeEnough));
  // TODO find the point (in the correct half of the graph) where the slope is the same
  // retrigger the ramp there
  // use a multi step process, each time refining further
  // correct half means: same half of the ramp we are on
  // half means: where f''(x) == 0 (slope of the slope, the steepest point)
  // in case of a simple sine shaper, it's 0.5, so don't worry for now


  // changeRate = ((gainStep/gainStep')-1)
  changeRate =
    // ((gainStep:max(smallest)/(gainStep':max(smallest)))-1)
    ((gainStep:max(ma.EPSILON)/(gainStep':max(ma.EPSILON)))-1)
    // / (fullDif:max(ma.EPSILON))
    / duration
    * (warpedSine(shape,rawRamp+rawRamp):max(ma.EPSILON)/(warpedSine(shape,rawRamp):max(ma.EPSILON)))
    * running;

  intervention =
    ((abs(changeRate)> (maxCR/ma.SR))
     * (1-dirChange)
    )
  ;
  // TODO: better value
  // maxCR = hslider("maxCR", 1000, 1, 6000, 1)*10;
  maxCR = 10000;
  compare(start,end,compSlope) =
    (
      select2(bigger , start , middle)
    , select2(bigger , middle , end)
    , compSlope
    )
  with {
    bigger = compSlope>slope(middle);
    slope(x) =
      shapeDif(shape,x,rampStep)
      *(1/(1-warpedSine(shape,x)));
    middle = (start+end)*.5;
  };
  // each compare halves the range,
  // so if we want the perfect newRamp value,
  // for a SR of 192k we need at least that many steps
  // 2^18 = 262144
  // so we need 18 compares in total
  // at 48k, 13 seems to little, 14 works
  // test with shape minimal, so 0.3 and duration = (3/16)^2
  // (lower shapes give jumps in the phase anyway)
  // with the above settings, too low nr of compares gives a stuck ramp
  newRamp =
    (start,end)
  , shapeDif(shape,rawRamp,rampStep')
    * ((rawDif'/rawDif)/(1-warpedSine(shape,rawRamp-rampStep)))
    :seq(i, 14, compare)
    : ((+:_*.5),!) // average start and end, throw away the rest

  with {
    start = 0;
    end = 1;
  };

  maxDerTable(shape) =
    // hslider("start", 0.2, 0, 1, 0.01)
    // maxDerivative(rampStep,shape)
    // maxDerivative(1/SIZE,shape)
    // ba.tabulate(0, maxDerivative(1/SIZE), SIZE, 0, 1, shape).val
    0
    // rdtable(SIZE,maxDerivative(1/SIZE,ba.time/SIZE),int(shape*SIZE))
    // + rampStep
    // 100/ma.SR
  with {
    SIZE = 1<<9;
    // SIZE = 1<<10 gives ocasional error values, presumably cause the dif becomes too small
    // more than 10 gives only errors
  };

  maxDerivative(stepsize,shape) =
    (0,1,shape,stepsize)
    : seq(i, 32, compareDer)//32 is overkill, but it's for a table, so it's OK
    : ((+:_*.5),!,!) // average start and end, throw away the rest
      * (shape!=0) // TODO: will this bite me in the *ss later on? is it even needed?
      // in any case, without this, it spits out a too high value.
    :select2(shape>(1-ma.EPSILON),_,1)
  ;
  compareDer(start,end,shape,stepsize) =
    (
      select2(bigger , start , middle)
    , select2(bigger , middle , end)
    , shape
    , stepsize
    )
  with {
    bigger = secondDerivative(shape,middle) > 0;
    derivative(shape,x) = warpedSine(shape,x+stepsize)-warpedSine(shape,x);
    secondDerivative(shape,x) = derivative(shape,x+stepsize)-derivative(shape,x);
    middle = (start+end)*.5;
  };


  compareShape(start,end,dif,compSlope) =
    (
      select2(bigger , start , middle)
    , select2(bigger , middle , end)
    , dif
    , compSlope
    )
  with {
    bigger = compSlope>slope;
    slope =
      (warpedSine(middle,x)-warpedSine(middle,x-rampStep))
      *(dif/(1-warpedSine(middle,x)));
    middle = (start+end)*.5;
    x = maxDerTable(middle);
  };
  inewShape =
    0.5;tmp=
          (start,end,rawDif)
        , (warpedSine(shape,rawRamp-rampStep)-warpedSine(shape,rawRamp-(2*rampStep)))
          * (rawDif'/(1-warpedSine(shape,rawRamp-rampStep)))
          :seq(i, 14, compareShape)
          : ((+:_*.5),!,!) // average start and end, throw away the rest
            // * ((ramp<(thres*2)) | ((ramp > (0.5-thres)) & (ramp < (0.5+thres))  ))
            // * releasing

            // <: select2(_<=0.0001,_,-1)
            // <: select2(_>=(1-0.0001),_,-1)
            // : select2(releasing,-1,_)
        with {
    start = 0.1;
    // end = 0.5;
    end = 0.9;
    // thres = (hslider("thres", 1, 0, 2, 0.01)/ma.SR)*10;
    // end = shape;
  };
  inewRampShapeX = maxDerTable(newShape');

  shapeIntervention =
    intervention
    &(ramp<thres)
    * releasing ;

  shapeInterventionHold =
    loop~_ with {
    loop(prev) =
      select2((shapeIntervention| (shapeIntervention'))
             , (prev*((ramp>ramp') ))
             , 1);
  };

  thres = hslider("thres", 1, 0, 1, 0.01)/rampStep/ma.SR*0.1;
  rampShapeFix =
    select2(shapeIntervention
           , ramp
           , newRampShapeX
           );
  switchedShape =
    select2(shapeInterventionHold
           , shapeSliderVal(shapeSlider)
           , ba.latch(shapeIntervention,newShape')
           );

  // };
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
  warpedSine(shape,x) =
    // the raw formula is faster than the tabulated version
    // 5 to 6 % CPU
    // warpedSineFormula(shape,x)
    par(i, nrShapes+1, table(i) * xfadeSelector(shapeSlider,i)):>_
  with {
    // 4.5 to 5.5 % CPU
    // 23 goes wrong with rel=3, shape=minimum, ramp not steep enough
    // SIZE = 1<<24;
    // table(i) = ba.tabulate(0, warpedSineFormula(shapeSliderVal(i)), SIZE, 0, 1, x).val;
    // 3 to 4 % CPU
    SIZE = 1<<14;
    table(i) = ba.tabulate(0, warpedSineFormula(shapeSliderVal(i)), SIZE, 0, 1, x).lin;
    // 4.5 to 6 % CPU
    // SIZE = 1<<9;
    // table(i) = ba.tabulate(0, warpedSineFormula(shapeSliderVal(i)), SIZE, 0, 1, x).cub;
    xfadeSelector(sel,nr) =
      ((sel<=nr)*((sel-nr)+1):max(0)) + ((sel>nr)*((nr-sel)+1):max(0));
  };


  warpedSineFormula(shape,x) =
    sineShaper(warp(shape,knee,x)):pow(power)
  with {
    power = (4*shape/3)+(1/3);
    knee = min(2*shape,2-(2*shape));
  };
  shapeSlider =
    // half;
    hslider("shape", 0, 0-half, half, 0.1) + half;
  shapeSliderVal(shapeSlider) =
    shapeSlider
    / nrShapes
    * range
    + start
    : hbargraph("shapeBG", 0.3, 0.7)
  with {
    range = 2* (.5-start);
    // lower shapes then 0.3 give jumps in the phase at low durations (d < (3/16:pow(2)))
    // also they give stuck ramps at nr of compares < 14
    // shapeSliderVal(shapeSlider) = hslider("shape", 0.5, 0.30, 0.70, 0.01);
    start = 0.3;
  };
  nrShapes = 8;
  half = nrShapes*.5;

  shape =
    select2 (attacking, shapeSliderVal(shapeSlider),1-shapeSliderVal(shapeSlider));

};
};


PMI_FBFFcompressor_N_chan(strength,thresh,att,rel,knee,prePost,link,FBFF,meter,N) =
  si.bus(N) <: si.bus(N*2):
  (
    ((ro.interleave(N,2) : par(i,N*2,abs) :par(i,N,it.interpolate_linear(FBFF)) : PMI_compression_gain_N_chan_db(strength*(1+((FBFF*-1)+1)),thresh,att,rel,knee,prePost,link,N)),si.bus(N))
    : (ro.interleave(N,2) : par(i,N,(meter: ba.db2linear)*_))
  )
  ~ si.bus(N);

PMI_compression_gain_N_chan_db(strength,thresh,att,rel,knee,prePost,link,1) =
  PMI_compression_gain_mono_db(strength,thresh,att,rel,knee,prePost);

PMI_compression_gain_N_chan_db(strength,thresh,att,rel,knee,prePost,link,N) =
  par(i,N,PMI_compression_gain_mono_db(strength,thresh,att,rel,knee,prePost))
  <: (si.bus(N),(ba.parallelMin(N) <: si.bus(N))) : ro.interleave(N,2) : par(i,N,(it.interpolate_linear(link)));

PMI_compression_gain_mono_db(strength,thresh,att,rel,knee,prePost) =
  // PMI(PMItime) : ba.bypass1(prePost,moog_AR(att,rel)) : ba.linear2db : gain_computer(strength,thresh,knee) : ba.bypass1((prePost!=1),moog_AR(rel,att))
  PMI(PMItime) : ba.linear2db : gain_computer(strength,thresh,knee)
                                // : ba.bypass1(prePost,ba.db2linear:my_AR(rel,att,r1,a1):ba.linear2db)
                                // : ba.bypass1((1-prePost),my_AR(rel,att,r1,a1))
  : ba.bypass1(prePost,ba.db2linear:my_AR(rel,att,r1,a1):ba.linear2db)
  : ba.bypass1((1-prePost),my_AR(rel,att,r1,a1))
    // : my_AR(0,att,0,a1)
    // : my_AR(rel,0,r1,0)
with {
  gain_computer(strength,thresh,knee,level) =
    select3((level>(thresh-(knee/2)))+(level>(thresh+(knee/2))),
            0,
            ((level-thresh+(knee/2)) : pow(2)/(2*max(ma.EPSILON,knee))),
            (level-thresh))
    : max(0)*-strength;
  PMI(time) =
    slidingPMI(s,192000,power)
  with {
    s = ba.sec2samp(time):int:max(1);
  };

};

my_AR(att,rel,a2,r2) = si.onePoleSwitching(att,rel) : der~_ with {
  // der(prev,x) = (x-prev) : (si.onePoleSwitching(a1,r1)+prev);
  der(prev,x) = (x-prev)*t(prev,x)+prev;
  t(prev,x) = select2(x>prev,a1,r1)*0.0001;
};

moog_AR(att,rel) = loop~_ with {
  loop(prev,x) = x:ve.moog_vcf_2bn(res,fr(x,prev)):ve.moog_vcf_2bn(res,fr(x,prev)):ve.moog_vcf_2bn(res,fr(x,prev)):ve.moog_vcf_2bn(res,fr(x,prev));
  fr(x,prev) = select2(x<=prev,1/att,1/rel);
  res = hslider("res", 0.1, 0.01, 10, 0.01);
};

slidingPMI(n,maxn,p) = pow(p) : ba.slidingMeanp(n,maxn) : pow(1/p);

AB(p) = ab:hgroup("[1]A/B",sel(aG(p),bG(p)));
sel(a,b,x) = select2(x,a,b);
aG(x) = vgroup("[0]a", x);
bG(x) = vgroup("[1]b", x);

N = 2;
B = si.bus(2);
ab = checkbox("[0]a/b");
bypass = AB(bypassP);
bypassP = checkbox("[00]bypass");
prePost = AB(prePostP);
prePostP = checkbox("[01]prePost");
strength = AB(strengthP);
strengthP = hslider("[02]strength", 20, 0, 100, 1) * 0.01;
thresh = AB(threshP);
threshP = hslider("[03]thresh",0,-30,6,1);
att = AB(attP);
attP = hslider("[04]attack",20,0,100,1)*0.001;
rel = AB(relP);
relP = hslider("[05]release",200,1,1000,1)*0.001;
knee = AB(kneeP);
kneeP = hslider("[06]knee",6,0,30,1);
link = AB(linkP);
linkP = hslider("[07]link", 60, 0, 100, 1) *0.01;
FBFF = AB(FBFFP);
FBFFP = hslider ("[08]fb-ff",50,0,100,1) *0.01;
power = AB(powerP);
powerP = hslider("[09]power", 2, ma.EPSILON, 10, 0.1);
PMItime = AB(PMItimeP);
PMItimeP = hslider("[10]PMI time",20,0,1000,1)*0.001;
dw = AB(dwP);
dwP = hslider ("[11]dry/wet",100,0,100,1) * 0.01:si.smoo;
a1 = AB(a1P);
a1P = hslider("[12]der attack",20,1,100,1);
r1 = AB(r1P);
r1P = hslider("[13]der release",200,1,1000,1);

ARtest = toggle(soft,loud) with {
  toggle(a,b) = select2(block,b,a);
  block = os.lf_sawpos(1)>0.5;
  soft = sine*0.001;
  loud = sine;
  sine = os.osc(5000);
};


meter =
  _<: attach(_, (max(-12):min(0):hbargraph(
                   "v:[10]meters/[unit:dB]", -12, 0)
                ));

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
