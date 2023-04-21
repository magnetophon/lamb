declare name "lamb";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPLv3";

import("stdfaust.lib");

// make a release based on:
// place

process =
  test:AR
       // test
       // PMI_FBFFcompressor_N_chan(strength,thresh,att,rel,knee,prePost,link,FBFF,meter,N);
       // (ARtest:PMI_compression_gain_mono_db(strength,thresh,att,rel,knee,prePost):ba.db2linear)
       // , os.lf_sawpos(1)>0.5
;


test0 = select3(
          ((os.lf_sawpos(1)>0.3)+(os.lf_sawpos(1)>0.5)),
          -0.1,0,1);
test1 = select2(os.lf_sawpos(1)>0.5, 0.1,0.9);
test =
  ((loop~_)
  , no.lfnoise(hslider("rate", 100, 0.1, 1000, 0.1))
  )
  :it.interpolate_linear(hslider("Xfade", 0, 0, 1, 0.1))
with {
  loop(prev,x) = no.lfnoise0(abs(prev*69)%9:pow(2)+1);
};


AR = loop~(_,_)
          // :(!,si.bus(4))
with {
  loop(prevRamp,prevGain,x) =
    ramp
   ,gain
   ,x
   , newShape
     // , newRampShapeX
     // ,(compareShape(0,1,0.5):>_)
     // ,abs(changeRate)
     // ,(changeRate* (warpedSine(shape,rawRamp+rampStep)/(warpedSine(shape,rawRamp):max(smallest))))
     // ,((maxCR/ma.SR))
     // , running
   , (maxDerTable(shape) :hbargraph("MD", 0, 1))
     // , maxDerivative(ba.time/(1<<16))
  with {
  rawRamp = (prevRamp+rampStep)*running:min(1):max(0);
  ramp = select2(intervention,rawRamp,newRamp)*running;
  // ramp = rawRamp;
  rampStep = 1 / ma.SR / duration;
  duration = select3(attacking+releasing*2,1,attack,release);
  attack = hslider("attack", 0.1, 0, 1, smallest):max(1 / ma.SR);
  release = hslider("release", 0.5, 0, 1, smallest):max(1 / ma.SR);
  // TODO better value
  smallest = 1/192000;
  gain = prevGain+gainStep:max(-1):min(1);
  rawGainStep = (shapedRamp-warpedSine(shape,rawRamp-rampStep))*fullDif;
  gainStep = select2(rawGainStep>0
                    , rawGainStep:min(0-smallest)
                    , rawGainStep:max(smallest)
                    )
             * running;
  rawDif = x-prevGain;
  fullDif =rawDif/(1-shapedRamp);
  running = (attacking | releasing) * (1-dirChange);
  dirChange = (attacking != attacking')| (releasing != releasing');
  // TODO: find proper N (needs to be bigger than 2 when compiling to 32 bit)
  // N = hslider("N", 1, 0, 8, 0.01);
  N=1;
  releasing = rawDif>(N / ma.SR);
  attacking = rawDif< 0-(N / ma.SR);
  // TODO find the point (in the correct half of the graph) where the slope is the same
  // retrigger the ramp there
  // use a multi step process, each time refining further
  // correct half means: same half of the ramp we are on
  // half means: where f''(x) == 0 (slope of the slope, the steepest point)
  // in case of a simple sine shaper, it's 0.5, so don't worry for now


  shapedRamp = warpedSine(shape,rawRamp);
  // changeRate = ((gainStep/gainStep')-1)
  changeRate = ((gainStep:max(smallest)/(gainStep':max(smallest)))-1)
               // / (fullDif:max(ma.EPSILON))
               / duration
               * (warpedSine(shape,rawRamp+rawRamp):max(ma.EPSILON)/(warpedSine(shape,rawRamp):max(ma.EPSILON)))
               * releasing;

  intervention =
    // 0;tmp=
    (abs(changeRate)> (maxCR/ma.SR))
    * (rawRamp > 0.01)
  ;
  // TODO: better value
  // maxCR = hslider("maxCR", 300, 1, 6000, 1)*10;
  maxCR = 10000;
  compare(start,end,dif,compSlope) =
    (
      select2(bigger , start , middle)
    , select2(bigger , middle , end)
    , dif
    , compSlope
    )
  with {
    bigger = compSlope>slope(middle);
    slope(x) =
      (warpedSine(shape,x)-warpedSine(shape,x-rampStep))
      *(dif/(1-warpedSine(shape,x)));
    middle = (start+end)*.5;
  };
  // each compare halves the range,
  // so if we want the perfect newRamp value,
  // for a SR of 192k we need at least that many steps
  // 2^18 = 262144
  // so we need 18 compares in total
  // at 48k, 13 total seems to little, 14 works
  newRamp =
    (start,end,rawDif)
  , (warpedSine(shape,rawRamp-rampStep)-warpedSine(shape,rawRamp-(2*rampStep)))
    * (rawDif'/(1-warpedSine(shape,rawRamp-rampStep)))
    :seq(i, 14, compare)
    : ((+:_*.5),!,!) // average start and end, throw away the rest

  with {
    start = 0;
    end = 0.5;
  };

  maxDerTable(shape) =
    maxDerivative(rampStep,shape)
    // ba.tabulate(0, maxDerivative(1/SIZE), SIZE, 0, 1, shape).val
    // rdtable(SIZE,maxDerivative(1/SIZE,ba.time/SIZE),int(shape*SIZE))
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


  compareShape(start,end,compSlope) =
    (
      select2(bigger , start , middle)
    , select2(bigger , middle , end)
    , compSlope
    )
  with {
    bigger = compSlope>=slope;
    slope =
      (warpedSine(middle,x+(1/ma.SR))-warpedSine(middle,x))
      *(1/(1-warpedSine(middle,x)));
    middle = (start+end)*.5;
    x = maxDerTable(middle);
  };
  newShape =
    (start,end)
  , (warpedSine(shape,rawRamp+(1/ma.SR))-warpedSine(shape,rawRamp))
    * (rawDif'/rawDif/(1-warpedSine(shape,rawRamp+(1/ma.SR))))
    :seq(i, 14, compareShape)
    : ((+:_*.5),!) // average start and end, throw away the rest
      // * ((ramp==0) | (ramp == 0.5))
      * running

      // <: select2(_<=0.0001,_,-1)
      // <: select2(_>=(1-0.0001),_,-1)
      // : select2(releasing,-1,_)
  with {
    start = 0.0001;
    end = 0.999;
  };
  newRampShapeX = maxDerTable(newShape);

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
  // sineShaper(x) = x:max(0.5):min(0.5);
  // sineShaper(x) = (0.5);
  warpedSine(shape,x) =
    sineShaper(warp(shape,knee,x)):pow(power)
  with {
    power = (4*shape/3)+(1/3);
    knee = min(2*shape,2-(2*shape));
  };
  shape = hslider("shape", 0.5, 0, 1, 0.01);
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
