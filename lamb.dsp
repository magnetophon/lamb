declare name "lamb";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPLv3";

import("/home/bart/source/lamb/stdfaust.lib");


process =
  // PMI_FBFFcompressor_N_chan(strength,thresh,att,rel,knee,prePost,link,FBFF,meter,2);
  //
  hgroup("",
         vgroup("[2]test", test)
         :vgroup("[1]AR",
                 AR
                ));

AR = loop~(_,_)
          // :(!,si.bus(4))
          // :(!,_)
with {
  loop(prevRamp,prevGain,x) =
    ramp
  , gain
  , x
  , (x:seq(i, 3, si.onePoleSwitching(releaseOP,attackOP)))
    // , (x==gain)
  with {
  duration =
    // select3(attacking+releasing*2,1,attack,release);
    (attack*attacking)+(release*releasing);
  attack = hslider("[1]attack time[scale:log]", 8, 0, 1000, 0.1)*0.001;
  release = hslider("[3]release time[scale:log]", 250, 0, 1000, 0.1)*0.001;
  attackOP = hslider("[5]OP attack time[scale:log]", 8, 0, 1000, 0.1)*0.001;
  releaseOP = hslider("[6]OPrelease time[scale:log]", 250, 0, 1000, 0.1)*0.001;
  gain = prevGain+gainStep ;
  gainStep =
    select2(releasing
           , rawGainStep :max(dif)
           , rawGainStep :min(dif)
           ) with {
    rawGainStep =
      shapeDif(shapeSlider,ramp,duration,ma.SR)*fullDif;
    fullDif =dif/(1-warpedSine(shapeSlider,ramp));
  };
  shapeDifFormula(shapeSlider,phase,len) =
    warpedSineFormula(shapeSlider,phase+len)
    - warpedSineFormula(shapeSlider,phase);

  // shapeDif(shape,phase,duration,sr) =
  // ba.tabulateNd(3,1,shapeDifFormula,nrShapes,1<<17,1<<7,0,0,1/192000/1,nrShapes,1,1/24000/(1/192000),shapeSlider,phase,(1 / ma.SR / duration));
  shapeDif(shapeSlider,phase,duration,sr) =
    warpedSine(shapeSlider,phase+(1 / sr / duration))
    - warpedSine(shapeSlider,phase);
  // warpedSineFormula(shapeSlider,phase+(1 / sr / duration))
  // - warpedSineFormula(shapeSlider,phase);

  dif = x-prevGain;
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
    :seq(i, 21, compare)
    : ((+:_*.5),!) // average start and end, throw away the rest
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
    ba.tabulateNd(1, warpedSineFormula,(nrShapes, SIZE,0, 0,nrShapes, 1, shapeSlider,x)).lin
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
    // table(i) = ba.tabulate(0, warpedSineFormula(shapeSliderVal(i)), SIZE, 0, 1, x).cub;
  };


  warpedSineFormula(shapeSlider,x) =
    sineShaper(warp(shape,knee,x)):pow(power)
  with {
    power = (4*shape/3)+(1/3);
    knee = min(2*shape,2-(2*shape));
    shape = shapeSliderVal(shapeSlider);
  };
  shapeSlider =
    // select2(releasing, 1-slider)
    select2(releasing
           , half+hslider("[2]attack shape" , 0, 0-half, half, 0.1)
           , half+hslider("[4]release shape", 0, 0-half, half, 0.1));

  nrShapes = 8;
  half = nrShapes*.5;

  shapeSliderVal(shapeSlider) =
    shapeSlider
    / nrShapes
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
  : ba.bypass1(prePost,ba.db2linear:AR:ba.linear2db)
  : ba.bypass1((1-prePost),AR)
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
  block = os.lf_sawpos(0.5)>0.5;
  soft = sine*0.1;
  loud = sine;
  sine = os.osc(5000);
};


meter =
  _<: attach(_, (max(-12):min(0):hbargraph(
                   "v:[10]meters/[unit:dB]", -12, 0)
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


mysel(x)=         (checkbox("AR")*
                   (si.onePoleSwitching(hslider("rel simple", 8, 0, 1000, 0.1)*0.001
                                       ,hslider("att simple", 8, 0, 1000, 0.1)*0.001,x))
                  , (1-checkbox("AR"))*AR(x):(!,_,!)):>_,x;
