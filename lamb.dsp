declare name "lamb";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPLv3";

import("stdfaust.lib");

///////////////////////////////////////////////////////////////////////////////
//                          compile time variables                           //
///////////////////////////////////////////////////////////////////////////////


// maxSampleRate = 192000;
maxSampleRate = 48000;
// maxSampleRate = 100; // for svg
// This has impact on the compile time and the cpu usage

nrChannels = 2;
// the number of input and output channels

enableGRout = 1;
// enable gain reduction outputs

selectSmoother = 0;
// 0 = just the sophisticated smoother, heavy on the CPU, long compile time
// 1 = just a regular 4-pole smoother with lookahead
// 2 = switchable between the two

selectConfiguration = 2;
// 0 = just the peak limiter
// 1 = just the serialGains
// 2 = both
// 3 = both, plus a leveler. TODO: FF/FB
// 4 = a "debugger" for the smoother algo, only usefull for developing
//
// the serialGains is "nrComps" serial compressors, ranging from slow to fast
// the serialGains GR is the sum of all of them
//
// when selectConfiguration = 2, the final GR is the minimum of the serialGains GR
// and the un-smoothed GR as determined by the limiter strength, threshold and gain
// the result of this is passed trough the lookahead smoother

nrComps = 8;
// the number of serial compressors

enableAB = 0;
// an A/B comparison system
// allows you to switch between two sets of parameters

enableDiffMeters = 0;
// a meter that shows you more or less GR of the peak limiter
// it can not show just the fast GR, since the smoother interacts with serialGains

serialGainsVarOrder = 1;
// the order of the smoothers in the serialGains is variable


maxOrder = 8;

///////////////////////////////////////////////////////////////////////////////
//                                  process                                     //
///////////////////////////////////////////////////////////////////////////////

process =
  lambSel(selectConfiguration);

lambSel(0) =
  (limiterGroup(selectConfiguration,preGain))
  : lamb(1)
  : postProc(enableGRout);

lambSel(1) =
  preGain<:
  (
    ( serialGainsGroup(selectConfiguration,serialGains)
      <: si.bus(nrChannels))
  , si.bus(nrChannels))
  : ro.interleave(nrChannels,2)
  : par(i, nrChannels, *);

lambSel(2) =
  preGain
  <:( (serialGainsGroup(selectConfiguration,serialGains)
      , si.bus(nrChannels))
      : lamb
    )
  : postProc(enableGRout);

lambSel(3) =
  preGain
  <:(
  leveler
  <:
  (serialGainsGroup(selectConfiguration,serialGains)
  , si.bus(nrChannels))
  : lamb)
  ~FBproc
   : postProc(enableGRout);

FBproc =
  ba.parallelMax(nrChannels)
, par(i, nrChannels, !)
;

leveler(prevGain) =
  ( si.bus(nrChannels)
    <: ( si.bus(nrChannels)
       , (levelerGain<:si.bus(nrChannels))))
  : ro.interleave(nrChannels,2)
  : par(i, nrChannels, *)
with {
  levelerGain =
    (avgGain - maxAvg)
    : max(0)
      *-1
    : levelerGroup(hbargraph("[99]leveler gain reduction", -24, 0))
    : ba.db2linear;
  att = levelerGroup(hslider("[01]leveler attack", 400, 1, 1000, 1))*0.001;
  rel = levelerGroup(hslider("[02]leveler release", 2000, 1, 10000, 1))*0.001;

  avgGain =
    // par(i, nrChannels, abs)
    // : ba.parallelMax(nrChannels)
    prevGain
    : si.onePoleSwitching(att,rel)
    : ba.slidingRMSp(rmsSamples,maxSampleRate)
    :ba.linear2db;
  rmsSamples = ba.sec2samp(rmsTime);
  rmsTime = levelerGroup(hslider("[03]leveler rms[unit:ms]", 400, minAttTimeMs, 1000, 1))*0.001;
  maxAvg = levelerGroup(hslider("[04]leveler threshold[unit:dB]", -3, -30, 30, 0.1));
};


lambSel(4) =
  SIN_tester;

lamb(parGain) =
  limiterGroup(selectConfiguration,lookahead_compressor_N_chan(parGain,strength,thresh,attack,release,knee,link,nrChannels));

preGain =
  par(i, nrChannels, _*inputGain);

///////////////////////////////////////////////////////////////////////////////
//                         SIN  smoother                                     //
///////////////////////////////////////////////////////////////////////////////

attackSamples = ba.sec2samp(attack);
maxAttackSamples =
  maxAttack*maxSampleRate
;

SIN(attack,release) = loop~(_,_)
with {
  loop(prevRamp,prevGain,x) =
    ramp
  , gain
  with {
  duration =
    (attack*attacking)+(release*releasing);
  gain = (prevGain+gainStep):min(x@attackSamples);
  gainStep =
    select2(releasing
           , rawGainStep :max(dif)
           , rawGainStep :min(dif)
           ) with {
    rawGainStep =
      shapeDif(shapeSlider,ramp,duration,ma.SR)*fullDif;
    fullDif =dif/(1-warpedSine(releasing,shapeSlider,ramp));
  };
  shapeDifFormula(shapeSlider,phase,len) =
    warpedSineFormula(shapeSlider,phase+len)
    - warpedSineFormula(shapeSlider,phase);

  shapeDif(shape,phase,duration,sr) =
    warpedSine(releasing,shapeSlider,phase+(1 / sr / duration))
    - warpedSine(releasing,shapeSlider,phase);

  hold = ba.slidingMin(attackSamples+1,maxAttackSamples,x);
  dif =
    hold
    - prevGain;

  releasing =
    dif>0;
  attacking =
    dif<0;

  ramp =
    (start,end)
  , shapeDif(shapeSlider,prevRamp+rampStep,duration',ma.SR)
    * ((dif'/dif)/(1-warpedSine(releasing,shapeSlider',prevRamp)))
    :seq(i, 16, compare)
    : ((+:_*.5),!) // average start and end, throw away the rest
    :max(0):min(1)
  with {
    start = 0;
    end = 1;
    rampStep = 1 / ma.SR / duration;

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
        *(1/(1-warpedSine(releasing,shapeSlider,x)));
      middle = (start+end)*.5;
    };
  };
  // ******************************************** the curves: ******************************


  warpedSine(releasing,shapeSlider,x) =
    newCurve(releasing,shapeSlider,x);
  // select2(checkbox("new")
  // select2(1
  // , OLDwarpedSine(releasing,shapeSlider,x)
  // , newCurve(releasing,shapeSlider,x)
  // );

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

  OLDwarpedSine(releasing,shapeSlider,x) =
    ba.tabulateNd(0, warpedSineFormula,(nrShapes, 1<<16,0, 0,nrShapes, maxRelease, shapeSlider,x)).lin;

  warpedSineFormula(shapeSlider,x) =
    sineShaper(warp(shape,knee,x:max(0):min(1))):pow(power)
  with {
    power = (4*shape/3)+(1/3);
    knee = min(2*shape,2-(2*shape));
    shape = shapeSliderVal(shapeSlider);
  };
  shapeSlider =
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
    start = 0.3;
  };
};
};

// *************************************** the NEW curves: ******************************

// Based on an algorithm by Dario Sanfilippo:
// https://www.desmos.com/calculator/2hxvf9q194
// Adapted by Bart Brouns:
// https://www.desmos.com/calculator/ubmqgogu2s
// simplified:
// https://www.desmos.com/calculator/cog4ujr7cs

f(x0,x1,y0,y1,k,x) = y0+(y1-y0)
                     * (1-exp((k*(x-x0) /(x1-x0))))
                     / (1-exp(k)) ;
fm1m2(c,x) =
  f(0,1,0,1,-2.42*(c),x);
s(p) = sin(p*ma.PI+1.5*ma.PI)*0.5+0.5;
c2(c,x) = s(fm1m2(c,x));
CurveFormula(c,x) =
  select2(c==0
          // , c2(c,x)
         , c2(c:pow(1+0.42*c),x)
         , s(x)
         );
Curve(c,x) =
  // CurveFormula(c,x);
  ba.tabulateNd(0, CurveFormula,(nrShapes, 1<<16,0, 0,1, 1, c,x)).lin;

newCurve(releasing,c,x)=
  select2(releasing
         , Curve(c,x *-1+1 )
           *-1+1
         , Curve(c,x)
         );



///////////////////////////////////////////////////////////////////////////////
//                                compressor                             //
///////////////////////////////////////////////////////////////////////////////
lookahead_compressor_N_chan(parGain,strength,thresh,att,rel,knee,link,N) =
  si.bus(N) <: si.bus(N*2):
  (
    par(i, N, _@(attackSamples+relHoldSamples))
   ,((par(i,N,abs) : lookahead_compression_gain_N_chan(parGain,strength,thresh,att,rel,knee,link,N))
     <: si.bus(N*2)
    )
  )
  :((ro.interleave(N,2)
     : par(i,N, *))
   , (si.bus(N)
      : diffMeters(enableDiffMeters,N,parGain)));

diffMeters(0,N,parGain) =
  si.bus(N);
diffMeters(1,N,parGain) =
  par(i, N,
      _<: attach(_,
                 (ba.linear2db -(parGain:ba.slidingMin(attackSamples+1,maxAttackSamples):ba.linear2db))
                 :min(0):vgroup("[99]fast gain reduction[unit:dB]",  hbargraph("%i", -12, 0))
                )
     );

lookahead_compression_gain_N_chan(parGain,strength,thresh,att,rel,knee,link,1) =
  lookahead_compression_gain_mono(parGain,strength,thresh,att,rel,knee);

lookahead_compression_gain_N_chan(parGain,strength,thresh,att,rel,knee,link,N) =
  par(i, N, ba.linear2db)
  <: (si.bus(N),(ba.parallelMax(N) <: si.bus(N)))
  : ro.interleave(N,2)
  : par(i,N,(it.interpolate_linear(link))
       )
  : par(i,N,lookahead_compression_gain_mono(parGain,strength,thresh,att,rel,knee)) ;

lookahead_compression_gain_mono(parGain,strength,thresh,att,rel,knee) =
  gain_computer(strength,thresh,knee)
  : ba.db2linear
  : min(parGain)
  : (releaseHold~_)
  : smootherSel(selectSmoother)
with {
  releaseHold(prevGain,rawGR) =
    max(
      min(prevGain,rawGR@relHoldSamples)
    , releaseLookahead(rawGR));
  releaseLookahead(rawGR) = rawGR:ba.slidingMin(relHoldSamples,maxSampleRate);
};


gain_computer(strength,thresh,knee,level) =
  select3((level>(thresh-(knee/2)))+(level>(thresh+(knee/2))),
          0,
          ((level-thresh+(knee/2)) : pow(2)/(2*max(ma.EPSILON,knee))),
          (level-thresh))
  : max(0)*-strength;

smootherSel(0) =
  SIN(attack,release) :(!,_);
smootherSel(1) =
  ba.slidingMin(attackSamples+1,maxAttackSamples)
  : smoother(4, release, attack ) ;
smootherSel(2) =
  _<: select2(SINsmoo
             , smootherSel(0)
             , smootherSel(1));

///////////////////////////////////////////////////////////////////////////////
//                                 SerialGains                                   //
///////////////////////////////////////////////////////////////////////////////

serialGains =
  level
  // : fakeFeedback
  : compArray
    // : slowSmoother
  : attachMeter(hbargraph("[99]serial gains GR", -12, 0))
with {
  level =
    ba.parallelMax(nrChannels):ba.linear2db;
  // TODO: clip the level at lim thresh?
  fakeFeedback = _;

  compArray(level) =
    (strengths,threshs,atts,rels,ordersAtt,ordersRel,knees)
    : (( 0,level)
      , ro.interleave(nrComps,7))
      // : par(i, nrComps, gainSerial(i,x))
      // : ba.parallelMin(nrComps)
    : seq(i, nrComps, gainSeq(i))
    : (ba.db2linear,!)
  ;

  strengths =
    LinArray(bottomStrength,topStrength,nrComps);
  threshs =
    LinArray(bottomThres,topThres,nrComps);
  atts =
    LogArray(bottomAtt,topAtt,nrComps);
  rels =
    LogArray(bottomRel,topRel,nrComps);
  ordersAtt =
    LinArray(bottomOrderAtt,topOrderAtt,nrComps)
    : par(i, nrComps, _+0.5:floor);
  ordersRel =
    LinArray(bottomOrderRel,topOrderRel,nrComps)
    : par(i, nrComps, _+0.5:floor);
  varOrders =
    LogArray(bottomOrder,topOrder,nrComps)
    : par(i, nrComps, _+0.5:floor);
  knees =
    LinArray(bottomKnee,topKnee,nrComps);


  gainSeq(i,prevGR,level,strength,thresh,att,rel,orderAtt,orderRel,knee) =
    (prevGR+GR)
  , (level+GR)
  , si.bus((nrComps-1-i)*7)
  with {
    GR =
      gain_computer(strength,thresh,knee,level)
      : ba.db2linear
      : smootherSelect(serialGainsVarOrder,orderRel,orderAtt,rel,att)

      : (ba.linear2db:hgroup("[98] gain reduction", vbargraph("%i", -12, 0))) ;
  };

  gainSerial(i,level,strength,thresh,att,rel,orderAtt,orderRel,knee) =
    gain_computer(strength,thresh,knee,level)
    : ba.db2linear
    : smootherSelect(serialGainsVarOrder,orderRel,orderAtt,rel,att)
      // use for nrComps < 10
    : attachMeter(hgroup("[98] gain reduction", vbargraph("%i", -12, 0)));
  // use for nrComps > 10    Cannot be empty because faust will make up a name containing the numner
  // : attachMeter(hgroup("[98] gain reduction", vbargraph("GR", -12, 0)));
  postSmoother = smoother(postOrder,postAtt,postRel);
  postOrder = 3;
};


smoother(order, att, rel, xx) =
  smootherOrder(order,order, att, rel, xx);

smootherOrder(maxOrder,order, att, rel, xx) =
  smootherARorder(maxOrder,order, order, att, rel, xx);

smootherARorder(maxOrder,orderAtt, orderRel, att, rel, xx) =
  xx : seq(i, maxOrder, loop(i) ~ _)
with {
  loop(i,fb, x) = coeff(i) * fb + (1.0 - coeff(i)) * x
  with {
  cutoffCorrection(order) = 1.0 / sqrt(pow(2.0, 1.0 / order) - 1.0);
  coeff(i) =
    ba.if(x > fb, attCoeff(i), relCoeff(i) );
  attCoeff(i) =
    exp(-TWOPIT * cutoffCorrection(orderAtt) / max(ma.EPSILON, att))
    * (i<orderAtt);
  relCoeff(i) =
    exp(-TWOPIT * cutoffCorrection(orderRel) / max(ma.EPSILON, rel))
    * (i<orderRel);
  TWOPIT = 2 * ma.PI * ma.T;
};
};

smootherSelect(0,orderAtt, orderRel,att,rel) =
  // smoother(4, att, rel);
  smootherARorder(4,4, 1, att, rel);
smootherSelect(1,orderAtt, orderRel,att,rel) =
  smootherARorder(maxOrder,orderAtt, orderRel, att, rel);


postProc(0) = meters(selectConfiguration):si.bus(nrChannels),par(i, nrChannels, !);
postProc(1) = meters(selectConfiguration):si.bus(nrChannels*2);

meters(0) = combineGroup(metersV);
meters(1) = metersH;
meters(2) = metersH;
meters(3) = meters(2);

metersH =
  ( si.bus(nrChannels)
  , par(i, nrChannels, (meterH(i))));
metersV =
  ( si.bus(nrChannels)
  , par(i, nrChannels, (meterV(i))));
///////////////////////////////////////////////////////////////////////////////
//                           parameter arrays                                //
///////////////////////////////////////////////////////////////////////////////
LinArray(bottom,top,0) =   0:! ;
LinArray(bottom,top,nrElements) =     par(i,nrElements,   ((top-bottom)*(i/(nrElements-1)))+bottom);

// make a log array of values, from bottom to top
LogArray(bottom,top,nrElements) =     par(i,nrElements,   pow((pow((top/bottom),1/(nrElements-1))),i)*bottom);

shapedArray(bottom,top,shape,0) =   0:! ;
shapedArray(bottom,top, shape ,nrElements) =
  par(i,nrElements,
      (i/(nrElements-1))
      // :shaper(shape)
      *(top-bottom)
      +bottom
     );
// with {
// https://www.desmos.com/calculator/pn4myus6x4
shaper(s,x) = (x-x*s)/(s-x*2*s+1);
// };
///////////////////////////////////////////////////////////////////////////////
//                                    GUI                                   //
///////////////////////////////////////////////////////////////////////////////
// selectConfiguration
// 0 = just the peak limiter
// 1 = just the serialGains
// 2 = both
// 3 = a "debugger" for the smoother algo, only usefull for developing

levelerGroup(x) = combineGroup(vgroup("[0]leveler", x));

serialGainsGroup(0,x) = vgroup("[0]", x);
serialGainsGroup(1,x) = serialGainsGroup(0,x);
serialGainsGroup(2,x) = combineGroup(vgroup("[1]serial gains", x));
serialGainsGroup(3,x) = serialGainsGroup(2,x);

limiterGroup(0,x) = combineGroup(vgroup("[0]", x));
limiterGroup(1,x) = limiterGroup(0,x);
limiterGroup(2,x) = combineGroup(vgroup("[2]peak limiter", x));
limiterGroup(3,x) =limiterGroup(2,x);


combineGroup(x) = hgroup("[02]", x);

meterH(i) =
  attachMeter(
    hbargraph(
      "v:[99]gain reduction/%i[unit:dB]", -12, 0)
  );
meterV(i) =
  attachMeter(
    vbargraph(
      "h:[99]gain reduction/%i[unit:dB]", -12, 0)
  );

attachMeter(b) =
  _<: attach(_, (ba.linear2db
                 // :max(-12):min(0)
                 : b));

AB(0,p) = p;
AB(1,p) = ab:hgroup("[2]",sel(aG(p),bG(p)));
sel(a,b,x) = select2(x,a,b);
aG(x) = vgroup("[0]a", x);
bG(x) = vgroup("[1]b", x);

SINsmoo =
  AB(enableAB,checkbox("SIN / 4-pole smoother"));

ab = checkbox("[1]a/b");
inputGain = AB(enableAB,hslider("[01]input gain", 0, -24, 24, 0.1)):ba.db2linear:si.smoo;
strength = AB(enableAB,strengthP);
strengthP = hslider("[02]strength", 100, 0, 100, 1) * 0.01;
thresh = AB(enableAB,threshP);
threshP = hslider("[03]thresh",-1,-30,0,0.1);
attack = AB(enableAB,attackP);
attackP = hslider("[04]attack[unit:ms] [scale:log]",9, 0, maxAttack*1000,0.1)*0.001;
attackShape = AB(enableAB,attackShapeP);
// attackShapeP = half+hslider("[05]attack shape" , 2, 0-half, half, 0.1);
attackShapeP = hslider("[05]attack shape" , 0, 0, 1, 0.01);
release = AB(enableAB,releaseP);
releaseP = hslider("[06]release[unit:ms] [scale:log]",60,1,maxRelease*1000,1)*0.001;
releaseShape = AB(enableAB,releaseShapeP);
// releaseShapeP = half+hslider("[07]release shape" , -3, 0-half, half, 0.1);
releaseShapeP = hslider("[07]release shape" , 0.5, 0, 1, 0.01);
relHoldSamples = AB(enableAB,relHoldSamplesP);
relHoldSamplesP = hslider("[08]release hold[unit:ms]", 50,
                          minAttTimeMs, maxAttack*1000, 1)*0.001:ba.sec2samp;
knee = AB(enableAB,kneeP);
kneeP = hslider("[09]knee",1,0,30,0.1);
link = AB(enableAB,linkP);
linkP = hslider("[10]link", 0, 0, 100, 1) *0.01;

//************************************** serialGains **********************************************************
postAtt = hslider("post attack[ms]", 0, 0, 2000, 1)*0.001;
postRel = hslider("post release[ms]", 2000, 0, 20000, 1)*0.001;

BTgroup(x) = hgroup("[0]", x);
bottomGroup(x) = BTgroup(vgroup("[0]", x));
topGroup(x) = BTgroup(vgroup("[1]", x));

bottomStrength = bottomGroup(hslider("[02]slow strength",100,0,100,1)*0.01);
bottomThres = bottomGroup(hslider("[03]slow thresh",3,-30,30,0.1));
bottomAtt = bottomGroup(hslider("[04]slow attack[unit:ms] [scale:log]",100, 10, 3000,10)*0.001);
bottomOrderAtt = bottomGroup(hslider("[05]slow attack order", 4, 1, maxOrder, 1));
bottomRel = bottomGroup(hslider("[06]slow release[unit:s] [scale:log]",7000,50,10000,50)*0.001);
bottomOrderRel = bottomGroup(hslider("[07]slow release order", 1, 1, maxOrder, 1));
bottomKnee = bottomGroup(hslider("[08]slow knee",3,0,30,0.1));

topStrength = topGroup(hslider("[02]fast strength",100,0,100,1)*0.01);
topThres = topGroup(hslider("[03]fast thresh",-1,-30,30,0.1));
topAtt = topGroup(hslider("[04]fast attack[unit:ms] [scale:log]",1, minAttTimeMs, 300,0.1)*0.001);
topOrderAtt = topGroup(hslider("[05]fast attack order", 4, 1, maxOrder, 1));
topRel = topGroup(hslider("[06]fast release[unit:s] [scale:log]",300,10,10000,10)*0.001);
topOrderRel = topGroup(hslider("[07]fast release order", 1, 1, maxOrder, 1));
topKnee = topGroup(hslider("[08]fast knee",0,0,30,0.1));

minAttTimeMs = 1000/maxSampleRate;

// *******************************   more constants *********************************************************

nrShapes = 8;
half = (nrShapes-1)*.5;

// 50 ms
maxAttack = 0.05;
// 0.5 sec
maxRelease = 0.5:max(maxAttack);

SINtest = toggle(soft,loud) with {
  toggle(a,b) = select2(block,b,a);
  block = os.lf_sawpos(0.5)>0.5;
  soft = sine*0.1;
  loud = sine;
  sine = os.osc(5000);
};



///////////////////////////////////////////////////////////////////////////////
//                                    test                                   //
///////////////////////////////////////////////////////////////////////////////


SIN_tester =
  hgroup("",
         vgroup("[2]test", test)
         <:vgroup("[1]SIN",
                  (
                    // ba.slidingMin(attackSamples+1,maxAttackSamples):
                    SIN(attack,release))
                  ,_@attackSamples
                   // ,ba.slidingMin(attackSamples+1,maxAttackSamples)
                   // ,(((ba.slidingMin(attackSamples+1,maxAttackSamples):smootherCascade(4, release, attack )),_@attackSamples):min)
                 ));
test =
  (select3(hslider("test", 2, 0, 2, 1)
          , test0
          , test1
          , test2
          )
  , no.lfnoise(hslider("rate", 100, 0.1, 20000, 0.1))
  )
  :it.interpolate_linear(hslider("Xfade", 0, 0, 1, 0.01));

test0 = select2(os.lf_sawpos(0.5)>0.5, -1,1);
test1 = select3(
          ((os.lf_sawpos(1)>hslider("POS1", 0.25, 0, 1 , 0.01))+(os.lf_sawpos(1)>hslider("POS2", 0.5, 0, 1 , 0.01))),
          1, -1, 0);
test2 =
  (loop~_)
with {
  loop(prev,x) = no.lfnoise0(abs(prev*69)%9:pow(0.75)*5+1);
};
