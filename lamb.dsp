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
// This has a lot of impact on the compile time and the cpu usage

nrChannels = 2;
// the number of input and output channels

selectOutputs = 0;
// 0 = just nrChannels audio outputs
// 1 = nrChannels audio + nrChannels gain reduction outputs
// 2 = for lamb-rs:
//  - nrChannels audio outputs
//  - a max level output,
//  downsampled to one value every samplesPerPixel samples
//  selectable between pre- and post- limiter
//  time aligned to the output audio
//  - a gain reduction output,
//  downsampled to one value every samplesPerPixel samples

selectSmoother = 0;
// 0 = just the sophisticated smoother, heavy on the CPU, long compile time
// 1 = just a regular 4-pole smoother with lookahead
// 2 = switchable between the two

selectConfiguration = 0;
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

levVarOrder = 0;
// the order of the smoothers in the serialGains is variable
//
serialGainsVarOrder = 0;
// the order of the smoothers in the serialGains is variable

// Display a bargraph to show the current latency.
// Mainly for https://github.com/magnetophon/lamb-rs to be able to report the latency to the host.
enableLatencyMeter = 0;


maxOrder = 4;

///////////////////////////////////////////////////////////////////////////////
//                                  process                                     //
///////////////////////////////////////////////////////////////////////////////

process =
  configSelector(selectConfiguration);
// DJcomp;

configSelector(0) =
si.bus(nrChannels)
<:
si.bus(nrChannels)
, (
  preGain
  : lamb(meters(selectConfiguration),1)
)
  : postProc(nrChannels,selectOutputs)
;

configSelector(1) =
  preGain<:
  (
    // ( serialGainsGroup(selectConfiguration,serialGains(1))
    ( DJcompression_gain_N_chan(DJstrength,DJthresh,DJattack,fastRelease,DJknee,1,nrChannels)
      <: si.bus(nrChannels))
  , si.bus(nrChannels))
  : ro.interleave(nrChannels,2)
  : par(i, nrChannels, *);

configSelector(2) =
  si.bus(nrChannels)
  <:
  si.bus(nrChannels)
,
  (
    preGain
    <:( (
        DJcompression_gain_N_chan(DJstrength,DJthresh,DJattack,fastRelease,DJknee,1,nrChannels)
      , si.bus(nrChannels))
        : lamb(meters(selectConfiguration))
      ))
  : postProc(nrChannels,selectOutputs)
;

configSelector(3) =

  si.bus(nrChannels)
  <:
  si.bus(nrChannels)
,
  (
    preGain
    <:
    (
      (
        (
          leveler
          <:
          (
            serialGainsGroup(selectConfiguration,serialGains)
          )
        )
      , si.bus(nrChannels)
      )
      : lamb(meters(selectConfiguration))
    )
    ~FBproc)
  : postProc(nrChannels,selectOutputs) ;

FBproc =
  par(i, nrChannels, !)
, (ba.parallelMin(nrChannels):ba.linear2db)
;

///////////////////////////////////////////////////////////////////////////////
//                                process                                     //
///////////////////////////////////////////////////////////////////////////////

// go down if the serialGains plus limiter are doing more work then the deadZone
// (fastGR*-1)> deadZone
// go up if leveler is down
// both with variable speed
// attTime = inf when prevGain

leveler(fastGR) =
  (si.bus(nrChannels))
  <: (level,si.bus(nrChannels))
  : (levelerGroup(
        loop(fastGR)~_
      )
    , si.bus(nrChannels))
with {
  level =
    par(i, nrChannels, abs)
    : ba.parallelMax(nrChannels):ba.linear2db;
  loop(fastGR,prevLevGR,levelX) =
    levGR
  with {
    levGR =
      // (prevLevGR+diff)
      newGR
      // (diff)
      // :min(0)

      // :max(hslider("maxL", 0, -12, 12, 0.1))
      // :min(hslider("minL", 0, -12, 12, 0.1))
      // :attachMeter(hgroup("", vbargraph("leveler GR[unit:dB]", -24, 0)));
      : hgroup("", hbargraph("leveler GR[unit:dB]", -24, 0));

    newGR =
      (
        (
          localDif
          // fastGR
          // : min(0)
          : smootherSelect(levVarOrder,levOrder,levOrder,levRel,levAtt)
            // : smootherOrder(maxOrder,levOrder
            //                 , levRel
            //                   /(upSpeedFactor:max(ma.EPSILON))
            //                   // *1@(4*maxAttackSamples)
            //                 , levAtt
            //                   /(downSpeedFactor:max(ma.EPSILON))
            //                   // *1@(4*maxAttackSamples)
            //                )
        )
        // - prevLevGR
        // -deadDown
      )
      // *1@(8*maxAttackSamples)
    ;
    levOrder= hslider("order", 1, 1, maxOrder, 1);
    levRel =
      hslider("[5]levRel", 2, 0, 10, 0.1)
      /(upSpeedFactor:max(ma.EPSILON));
    levAtt =
      hslider("[1]levAtt", 1, 0, 5, 0.05)
      /(downSpeedFactor:max(ma.EPSILON));

    diff = select2(levReleasing
                  , down
                  , up)
    ;
    levReleasing =
      prevLevGR>fastGR;
    // (localDif*-1)<deadDown;
    down =
      prevLevGR+
      (
        (localDif-deadDown)
        :min(0)
         // :max(-1)
         * (downSpeed/ma.SR)
         * downSpeedFactor
      )
    ;
    localDif = fastGR-prevLevGR;

    up =
      // max(fastGR,prevLevGR)
      // prevLevGR+localDif
      // (prevLevGR)
      // +
      (localDif+deadUp)
      // : min(0)
      // :min(down)
      :max(prevLevGR)
       // :max(down)
       // : min(0)
       // Order(maxOrder,hslider("order", 1, 1, maxO, step),
      : smootherOrder(maxOrder,hslider("order", 1, 1, maxOrder, 1)
                      , hslider("upRel", 0, 0, 5, 0.001)
                        /(upSpeedFactor:max(ma.EPSILON))
                        *1@(4*maxAttackSamples)
                      , 0
                     )
    ;
    OLDup =
      prevLevGR+
      (upSpeed/ma.SR
       * upSpeedFactor) ;

    upSpeedFactor =
      (levelX
       +prevLevGR
       // -threshP
       +deadUp
      )
      :min(0)
       // :hbargraph("sum[unit:dB]", -90, 0)
       /((levKneeUp-deadUp):max(ma.EPSILON))
       *-1
      :min(1)
       * levReleasing
       // : si.onePoleSwitching(hslider("Art", 0, 0, 0.5, 0.001),0)
      : smoother(1,hslider("[7]Rellll", 60, 0, 500, 1)*0.001
                   // , hslider("Art", 0, 0, 0.5, 0.001)
                 , 0
                )
        // : Curve(hslider("curve", 0, 0, 1, 0.01))
      : hbargraph("[8]usf", 0, 1);

    downSpeedFactor = (localDif*-1-deadDown)/levKneeDown:max(0):min(1)
                                                                // *(1-levReleasing)
                      :hbargraph("[4]dsf", 0, 1);
    downSpeed = hslider("[1]down speed[unit:dB/S]", 2, 0, 6, 0.01);
    deadDown = hslider("[2]dead down[unit:dB]", 6, ma.EPSILON, 12, 0.1);
    levKneeDown = hslider("[3]levKnee down[unit:dB]", 24, ma.EPSILON, 30, 0.1);
    upSpeed = hslider("[5]up speed[unit:dB/S]", 2, 0, 6, 0.01);
    deadUp = hslider("[6]dead up[unit:dB]", 6, 0, 30, 0.1);
    levKneeUp = hslider("[7]levKnee Up[unit:dB]", 24, 0, 90, 0.1);
  };

  // hslider("lev gain", 0, -12, 12, 0.1),
  // si.bus(nrChannels);
};


//*********************************************************************************
OLDLeveler(prevGain) =
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
    : levelerGroup(hbargraph("[99]leveler gain reduction[unit:dB]", -24, 0))
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


configSelector(4) =
  SIN_tester;

lamb(meters,parGain) =
  (parGain,strength,thresh,attack,attackShape,release,releaseShape,knee,link,relHoldSamples,si.bus(nrChannels))
  :lookahead_compressor_N_chan(nrChannels,meters)
;

preGain =
  par(i, nrChannels, _*inputGain(selectConfiguration));

///////////////////////////////////////////////////////////////////////////////
//                         SIN  smoother                                     //
///////////////////////////////////////////////////////////////////////////////

maxAttackSamples =
  maxAttack*maxSampleRate;
maxRelHoldSamples =
  maxRelHold*maxSampleRate;
attackSamples = ba.sec2samp(attack);
fullLatency =
  (attackSamples+relHoldSamples)
  * (1-(selectConfiguration==1))
  :latencyMeter(enableLatencyMeter);
latencyMeter(0) = _;
latencyMeter(1) = hbargraph("[99]latency[unit:samples]", 0, maxAttackSamples+maxRelHoldSamples);

SIN(attack,attackShape,release,releaseShape) = loop~(_,_)
with {
  loop(prevRamp,prevGain,x) =
    ramp
  , gain
  with {
  duration =
    (attack*attacking)+(release*releasing);
  gain =
    (prevGain+gainStep)
    <:select2(attacking
             , max(_,prevGain)
             , min(_,prevGain))
    :min(x@attackSamples);

  gainStep =
    rawGainStep
    // select2(releasing
    // , rawGainStep :max(dif)
    // , rawGainStep :min(dif)
    // )
  with {
    rawGainStep =
      select2(duration==0
             , shapeDif(shapeSlider,ramp,duration,ma.SR)*fullDif
             , dif
             );
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
     // :seq(i, 24, compare)
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

  shapeSlider =
    select2(releasing
           , attackShape
           , releaseShape);
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
  ba.tabulateNd(1, warpedSineFormula,(nrShapes, 1<<16,0, 0,nrShapes, maxRelease, shapeSlider,x)).lin;

warpedSineFormula(shapeSlider,x) =
  sineShaper(warp(shape,knee,x:max(0):min(1))):pow(power)
with {
  power = (4*shape/3)+(1/3);
  knee = min(2*shape,2-(2*shape));
  shape = shapeSliderVal(shapeSlider);
};


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
         )
  :max(0):min(1)
;
Curve(c,x) =
  // CurveFormula(c,x);
  ba.tabulateNd(1, CurveFormula,(nrShapes, 1<<16,0, 0,1, 1, c,x)).lin;

newCurve(releasing,c,x) =
  select2(releasing
         , Curve(c,x *-1+1 )
           *-1+1
         , Curve(c,x)
         );



///////////////////////////////////////////////////////////////////////////////
//                                compressor                             //
///////////////////////////////////////////////////////////////////////////////
lookahead_compressor_N_chan(N,meters,parGain,strength,thresh,attack,attackShape,release,releaseShape,knee,link,relHoldSamples) =
  si.bus(N) <: si.bus(N*2):
  (
    par(i, N, _@fullLatency)
   ,((par(i,N,abs) : lookahead_compression_gain_N_chan(parGain,strength,thresh,attack,attackShape,release,releaseShape,knee,link,N,relHoldSamples))
     <: si.bus(N*2)
    )
  )
  :(
    (si.bus(N), meters)
    :(ro.interleave(N,2)
      : par(i,N, *))
  , (si.bus(N)
    )) ;

lookahead_compression_gain_N_chan(parGain,strength,thresh,attack,attackShape,release,releaseShape,knee,link,1,relHoldSamples) =
  lookahead_compression_gain_mono(parGain,strength,thresh,attack,attackShape,release,releaseShape,knee,relHoldSamples);

lookahead_compression_gain_N_chan(parGain,strength,thresh,attack,attackShape,release,releaseShape,knee,link,N,relHoldSamples) =
  par(i, N, ba.linear2db)
  <: (si.bus(N),(ba.parallelMax(N) <: si.bus(N)))
  : ro.interleave(N,2)
  : par(i,N,(it.interpolate_linear(link))
       )
  : par(i,N,lookahead_compression_gain_mono(parGain,strength,thresh,attack,attackShape,release,releaseShape,knee,relHoldSamples)) ;

lookahead_compression_gain_mono(parGain,strength,thresh,attack,attackShape,release,releaseShape,knee,relHoldSamples) =
  gain_computer(strength,thresh,knee)
  : ba.db2linear
  : min(parGain) // TODO: leave pargain in dB and do the min before the conversion here?
  : (releaseHold~_)
  : smootherSel(selectSmoother,attack,attackShape,release,releaseShape)
with {
  releaseHold(prevGain,rawGR) =
    max(
      min(prevGain,rawGR@relHoldSamples)
    , releaseLookahead(rawGR));
  releaseLookahead(rawGR) = rawGR:ba.slidingMin(relHoldSamples+1,maxRelHoldSamples);
};


gain_computer(strength,thresh,knee,level) =
  select3((level>(thresh-(knee/2)))+(level>(thresh+(knee/2))),
          0,
          ((level-thresh+(knee/2)) : pow(2)/(2*max(ma.EPSILON,knee))),
          (level-thresh))
  : max(0)*-strength;

smootherSel(0,attack,attackShape,release,releaseShape) =
  SIN(attack,attackShape,release,releaseShape) :(!,_);
smootherSel(1,attack,attackShape,release,releaseShape) =
  ba.slidingMin(attackSamples+1,maxAttackSamples)
  : smoother(4, release, attack )
with {
  // attackSamples = ba.sec2samp(attack);
};
smootherSel(2,attack,attackShape,release,releaseShape) =
  _<: select2(SINsmoo
             , smootherSel(0,attack,attackShape,release,releaseShape)
             , smootherSel(1,attack,attackShape,release,releaseShape));

///////////////////////////////////////////////////////////////////////////////
//                                 SerialGains                                   //
///////////////////////////////////////////////////////////////////////////////

serialGains(levelerGain) =
  level+levelerGain
  // : fakeFeedback
  : compArray
    // : slowSmoother
  : hbargraph("[99]serial gains GR[unit:dB]", -6, 0)
    + levelerGain
  :ba.db2linear
with {
  level =
    par(i, nrChannels, abs)
    : ba.parallelMax(nrChannels):ba.linear2db;
  // TODO: clip the level at lim thresh?
  fakeFeedback = _;

  compArray(level) =
    (strengths,threshs,atts,rels,ordersAtt,ordersRel,knees)
    : (( 0,level)
      , ro.interleave(nrComps,7))
      // : par(i, nrComps, gainSerial(i,x))
      // : ba.parallelMin(nrComps)
    : seq(i, nrComps, gainSeq(i))
    : (_,!)
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

      : (ba.linear2db:hgroup("[98] gain reduction", vbargraph("%i[unit:dB]", -6, 0))) ;
  };

  gainSerial(i,level,strength,thresh,att,rel,orderAtt,orderRel,knee) =
    gain_computer(strength,thresh,knee,level)
    : ba.db2linear
    : smootherSelect(serialGainsVarOrder,orderRel,orderAtt,rel,att)
      // use for nrComps < 10
    : attachMeter(hgroup("[98] gain reduction", vbargraph("%i[unit:dB]", -6, 0)));
  // use for nrComps > 10    Cannot be empty because faust will make up a name containing the numner
  // : attachMeter(hgroup("[98] gain reduction", vbargraph("GR", -12, 0)));
  postSmoother = smoother(postOrder,postAtt,postRel);
  postOrder = 3;
};

///////////////////////////////////////////////////////////////////////////////
//                               smoothers                                   //
///////////////////////////////////////////////////////////////////////////////

// smoother adapted from Dario Sanfilippo
// https://github.com/dariosanfilippo/limiterStereo/blob/da1c38cc393f08b5dd79e56ffd4e6256af07a708/limiterStereo.dsp#L90-L101

// fixed order
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
  smootherARorder(1,1, 1, att, rel);
smootherSelect(1,orderAtt, orderRel,att,rel) =
  smootherARorder(maxOrder,orderAtt, orderRel, att, rel);


// postProc(0) = meters(selectConfiguration):si.bus(nrChannels),par(i, nrChannels, _<:attach(!,_));
// postProc(0) = si.bus(nrChannels),par(i, nrChannels, !);
// postProc(1) = si.bus(nrChannels*2);

meters(0) = combineGroup(metersV);
meters(1) = metersH;
meters(2) = meters(1);
meters(3) = meters(1);

metersH =
  ( par(i, nrChannels, (meterH(i)))
  , si.bus(nrChannels)
  );
metersV =
  ( par(i, nrChannels, (meterV(i)))
  , si.bus(nrChannels)
  );
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
      "v:[99]gain reduction/%i[unit:dB]", -24, 0)
  );
meterV(i) =
  attachMeter(
    vbargraph(
      "h:[99]gain reduction/%i[unit:dB]", -24, 0)
  );

attachMeter(b) =
  _<: attach(_, (ba.linear2db
                 : b));

AB(0,p) = p;
AB(1,p) = ab:hgroup("[2]",sel(aG(p),bG(p)));
sel(a,b,x) = select2(x,a,b);
aG(x) = vgroup("[0]a", x);
bG(x) = vgroup("[1]b", x);

SINsmoo =
  AB(enableAB,checkbox("SIN / 4-pole smoother"));

ab = checkbox("[1]a/b");

inputGain(0) =
  limiterGroup(selectConfiguration,
               gainAB) ;
inputGain(selectConfiguration) =
  gainAB;

gainAB = AB(enableAB,hslider("[01]input gain", 0, -24, 24, 0.1)):ba.db2linear:si.smoo;
strength =
  limiterGroup(selectConfiguration,
               AB(enableAB,strengthP));
strengthP = hslider("[02]strength", 100, 0, 100, 1) * 0.01;
thresh =
  limiterGroup(selectConfiguration,
               AB(enableAB,threshP));
threshP = hslider("[03]thresh",-1,-30,0,0.1);
attack =
  limiterGroup(selectConfiguration,
               AB(enableAB,attackP));
attackP = hslider("[04]attack[unit:ms] [scale:log]",9, 0, maxAttack*1000,0.1)*0.001;
attackShape =
  limiterGroup(selectConfiguration,
               AB(enableAB,attackShapeP));
// attackShapeP = half+hslider("[05]attack shape" , 2, 0-half, half, 0.1);
attackShapeP = hslider("[05]attack shape" , 0, 0, 1, 0.01);
release =
  limiterGroup(selectConfiguration,
               AB(enableAB,releaseP));
releaseP = hslider("[06]release[unit:ms] [scale:log]",60,1,maxRelease*1000,1)*0.001;
releaseShape =
  limiterGroup(selectConfiguration,
               AB(enableAB,releaseShapeP));
// releaseShapeP = half+hslider("[07]release shape" , -3, 0-half, half, 0.1);
releaseShapeP = hslider("[07]release shape" , 0.5, 0, 1, 0.01);
relHoldSamples =
  limiterGroup(selectConfiguration,
               AB(enableAB,relHoldSamplesP));
relHoldSamplesP = hslider("[08]release hold[unit:ms]", 50,
                          minAttTimeMs, maxRelHold*1000, 1)*0.001:ba.sec2samp;
knee =
  limiterGroup(selectConfiguration,
               AB(enableAB,kneeP));
kneeP = hslider("[09]knee",1,0,30,0.1);
link =
  limiterGroup(selectConfiguration,
               AB(enableAB,linkP));
linkP = hslider("[10]link", 0, 0, 100, 1) *0.01;

//************************************** serialGains **********************************************************
postAtt = hslider("post attack[ms]", 0, 0, 2000, 1)*0.001;
postRel = hslider("post release[ms]", 2000, 0, 20000, 1)*0.001;

BTgroup(x) = hgroup("[0]", x);
bottomGroup(x) = BTgroup(vgroup("[0]", x));
topGroup(x) = BTgroup(vgroup("[1]", x));

bottomStrength = bottomGroup(hslider("[02]slow strength",100,0,100,1)*0.01);
bottomThres = bottomGroup(hslider("[03]slow thresh",3,-30,30,0.1));
bottomAtt = bottomGroup(hslider("[04]slow attack[unit:ms] [scale:log]",700, 10, 3000,10)*0.001);
bottomOrderAtt = bottomGroup(hslider("[05]slow attack order", 4, 1, maxOrder, 1));
bottomRel = bottomGroup(hslider("[06]slow release[unit:s] [scale:log]",700,50,5000,50)*0.001);
bottomOrderRel = bottomGroup(hslider("[07]slow release order", 1, 1, maxOrder, 1));
bottomKnee = bottomGroup(hslider("[08]slow knee",4,0,30,0.1));

topStrength = topGroup(hslider("[02]fast strength",100,0,100,1)*0.01);
topThres = topGroup(hslider("[03]fast thresh",-1,-30,30,0.1));
topAtt = topGroup(hslider("[04]fast attack[unit:ms] [scale:log]",9, minAttTimeMs, 300,0.1)*0.001);
topOrderAtt = topGroup(hslider("[05]fast attack order", 4, 1, maxOrder, 1));
topRel = topGroup(hslider("[06]fast release[unit:s] [scale:log]",60,10,1000,10)*0.001);
topOrderRel = topGroup(hslider("[07]fast release order", 1, 1, maxOrder, 1));
topKnee = topGroup(hslider("[08]fast knee",0,0,30,0.1));

minAttTimeMs = 1000/maxSampleRate;

// *******************************   more constants *********************************************************

nrShapes = 3;
half = (nrShapes-1)*.5;

// 50 ms
maxAttack = 0.05;
// 1 s
// maxAttack = 1;
// 0.5 sec
maxRelease = 0.5:max(maxAttack);
// 50 ms
maxRelHold = 0.05;

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

test0 = select2(os.lf_sawpos(hslider("block freq", 0.5, 0.1, 4, 0.1))>0.5, -1,1);
test1 = select3(
          ((os.lf_sawpos(hslider("block freq", 0.5, 0.1, 4, 0.1))>hslider("POS1", 0.25, 0, 1 , 0.01))+(os.lf_sawpos(1)>hslider("POS2", 0.5, 0, 1 , 0.01))),
          1, -1, 0);
test2 =
  (loop~_)
with {
  loop(prev,x) = no.lfnoise0(abs(prev*69)%9:pow(0.75)*5+1);
};

///////////////////////////////////////////////////////////////////////////////
//                                  notes                                    //
///////////////////////////////////////////////////////////////////////////////
// TODO: https://www.desmos.com/calculator/ubmqgogu2s
// old scaled: https://www.desmos.com/calculator/fa6rreh99a
// resample docs: https://henquist.github.io/2.0.0/

///////////////////////////////////////////////////////////////////////////////
//                                  anti-pump                                    //
///////////////////////////////////////////////////////////////////////////////

DJcomp =
  DJcompressor_N_chan(DJstrength,DJthresh,DJattack,fastRelease,DJknee,1,2) ;

DJcompressor_N_chan(strength,thresh,att,rel,knee,link,N) =
  par(i, N, _*inputGain(selectConfiguration))
  <: si.bus(N*2)
  : (
  si.bus(N)
 ,(DJcompression_gain_N_chan(strength,thresh,att,rel,knee,link,N)
   <: si.bus(N*2)
  )
)
  : (
    ro.interleave(N,2)
    : par(i,N, *)
  )
, si.bus(N)
  : postProc(N,selectOutputs)
  ;

  DJcompression_gain_N_chan(strength,thresh,att,rel,knee,link,1) =
    abs:ba.linear2db
    : DJcompression_gain_mono(strength,thresh,att,rel,knee);

  DJcompression_gain_N_chan(strength,thresh,att,rel,knee,1,N) =
    par(i, N, abs)
    : ba.parallelMax(N)
    : ba.linear2db
    : DJcompression_gain_mono(strength,thresh,att,rel,knee);

  DJcompression_gain_N_chan(strength,thresh,att,rel,knee,link,N) =
    par(i, N, abs:ba.linear2db)
    <: (si.bus(N),(ba.parallelMax(N) <: si.bus(N)))
    : ro.interleave(N,2)
    : par(i,N,(it.interpolate_linear(link))
         )
    : par(i,N,DJcompression_gain_mono(strength,thresh,att,rel,knee)) ;


  DJcompression_gain_mono(strength,thresh,att,rel,knee,level) =
    loop~(_,_)
         : (_,!)
         : ba.db2linear
         : smootherARorder(maxOrder, orderRel,orderAtt, 0, att)
  with {
    loop(prevGain,prevRef) =
      gain,ref
    with {
    gain =
      (  gain_computer(1,thresh,knee,level)
         : ba.db2linear
         : smootherARorder(maxOrder, orderRel,orderAtt, adaptiveRel, 0)
           // , ((level-limThres):max(0)*-1: ba.db2linear)
      )
      // :min(gain_computer(1,thresh+limThres,limKnee,level): ba.db2linear)
      // : smootherARorder(maxOrder, orderRelLim,4, releaseLim, 0)
      : ba.linear2db
        * strength
      : attachLatency(hbargraph("slow GR[unit:dB]", -24, 0))
    ;

    adaptiveRel =
      fade_to_inf(
        // shaper(adaptShape,
        1-dv
        // )
       ,rel) ;

    ref =
      // (prevGain+transitionRange)
      (prevGain-dvBot)
      // : min(transitionRange)
      : ba.db2linear
        // : smootherOrder(maxOrder,refOrder,refRel,0)
      : smootherOrder(1,1,refRel,0)
      : ba.linear2db
      : attachLatency(hbargraph("ref[unit:dB]", -24, 0))
    ;
    refRel =
      interpolate_logarithmic(
        // shaper(refShape,
        refDv
        // )
      , slowRelease,slowRelease/ma.EPSILON) ;
    refDv =
      it.remap(refBot, refTop, 1, 0,fastGR:min(refTop):max(refBot))
      : attachLatency(hbargraph("ref dv", 0, 1));
    dv =
      it.remap(dvBot, dvTop, 1, 0,fastGR:min(dvTop):max(dvBot))
      : attachLatency(hbargraph("dv", 0, 1));
    fastGR =
      (prevGain-prevRef)
      // :min(0)
      // :hbargraph("fast GR[unit:dB]", -24, 0)
      // :hbargraph("fast GR plus[unit:dB]", 0, 12 )
    ;
    attachLatency(m) = _<:attach(_,(_@fullLatency):m);
  };
  };


  // gain_computer(strength,thresh,knee,level) =
  // select3((level>(thresh-(knee/2)))+(level>(thresh+(knee/2))),
  // 0,
  // ((level-thresh+(knee/2)) : pow(2)/(2*max(ma.EPSILON,knee))),
  // (level-thresh))
  // : max(0)*-strength;

  postProc(N,0) =
    par(i, N, !),si.bus(N),par(i, N, !);
  postProc(N,1) =
    par(i, N, !),si.bus(N*2);
  postProc(N,2) =
    (
      si.bus(N)
    , ( si.bus(N)
        <:  si.bus(N*2)
      )
    , si.bus(N)
    )
    :(
    ro.crossnn(N)
  , si.bus(N*2)
  )
    : (
      si.bus(N)
    , par(i, 2, soloOp(max))
    , soloOp(min))
    : (
      si.bus(N)
    , (select2(inputSel
              ,_@fullLatency
              ,_)
       :abs:rm.reduce(max,samplesPerPixel))
    , rm.reduce(min,samplesPerPixel)
    )
  with {
    soloOp(op) =
      si.bus(N)
      <:(ba.parallelOp(op,N)
        , si.bus(N)
        )
      : ba.selectn(N+1,channelSelect);
    // 0 = all channels combined, 1 is first channel, 2 = second channel .. N = Nth channel
    channelSelect = hslider("[0]channel select", 0, 0, N, 1);
    inputSel = hslider("[1]input, output", 0, 0, 1, 1);
    samplesPerPixel = hslider("[2]samples per pixel", 16, 1, 1024, 1);
  };


  ///////////////////////////////////////////////////////////////////////////////
  //                              Utilities                                    //
  ///////////////////////////////////////////////////////////////////////////////

  fade_to_inf(dv,v0) =
    v0/max(1-dv,ma.EPSILON);

  interpolate_logarithmic(dv,v0,v1) =
    pow(((v1/v0)),dv)*v0 ;

  // https://www.desmos.com/calculator/pn4myus6x4
  shaper(s,x) = (x-x*s)/(s-x*2*s+1);
  ///////////////////////////////////////////////////////////////////////////////
  //                                    GUI                                   //
  ///////////////////////////////////////////////////////////////////////////////

  oneKnob = hslider("anti pump", 1, 0, 1, 0.01):si.smoo;

  DJinputGain =
    0;
  // (it.remap(0, 0.5, -9, 0,oneKnob:min(0.5))):ba.db2linear;
  DJstrength =
    // 1;
    // it.remap(0, 0.5, 0, 1,oneKnob:min(0.5));
    serialGainsGroup(selectConfiguration,
                     hslider("[02]DJ strength[unit:%]", 100, 0, 100, 1) * 0.01);

  DJthresh =
    thresh;
  // -1;
  DJattack =
    // 0;
    // 0.009;
    serialGainsGroup(selectConfiguration,
                     hslider("[04]DJ attack[unit:ms] [scale:log]",1, 1000/maxSampleRate, 300,0.1)*0.001);
  orderAtt =
    4;
  // hslider("[05]attack order", 4, 1, maxOrder, 1);
  fastRelease =
    release;
  // hslider("[06]fast release[unit:ms] [scale:log]",420,0.1,maxRelease*1000,1)*0.001;
  transitionRange =
    // it.remap(0, 1, 12, 6,oneKnob);
    // 9;
    serialGainsGroup(selectConfiguration,
                     hslider("[07]release transition range[unit:dB]",9,0,30,0.1));
  refTop =
    serialGainsGroup(selectConfiguration,
                     hslider("[07]ref top[unit:dB]",3,0,12,0.1));
  refBot =
    serialGainsGroup(selectConfiguration,
                     hslider("[07]ref bot[unit:dB]",-24,-30,0,0.1));
  dvTop =
    serialGainsGroup(selectConfiguration,
                     hslider("[07]dv top[unit:dB]",2,-6,6,0.1));
  dvBot =
    serialGainsGroup(selectConfiguration,
                     hslider("[07]dv bot[unit:dB]",-7,-30,0,0.1));
  slowRelease =
    serialGainsGroup(selectConfiguration,
                     hslider("[08]slow release[unit:ms] [scale:log]",1000,50,10000,50)*0.001);
  // it.remap(0, 1, 0.5, 4,oneKnob);
  orderRel =
    1;
  // hslider("[09]release order", 1, 1, maxOrder, 1);
  DJknee =
    knee;
  // hslider("[10]knee[unit:dB]",1,0,72,0.1);
  // it.remap(0, 1, 9, 9,oneKnob:max(0.5));
  refOrder =
    serialGainsGroup(selectConfiguration,
                     hslider("[11]ref release order", 1, 1, maxOrder, 1));

  // give it some headroom.
  // when built into a mixer, there should be a limiter on the master.
  // if that is done, this can come back up to -1dB
  limThres =
    serialGainsGroup(selectConfiguration,
                     hslider("[14]AP lim offset[unit:dB]",30,-30,30,0.1));
  releaseLim =
    // hslider("[15]release limiter[unit:ms] [scale:log]",60,5,500,1)*0.001;
    it.remap(0, 1, 0.12, 0.06,oneKnob);

  orderRelLim =
    4;
  // hslider("[16]lim release order", 4, 1, maxOrder, 1);
  limKnee =
    knee;
  // hslider("[17]lim knee[unit:dB]",1,0,72,0.1);
  // it.remap(0.5, 1, 12, 0,oneKnob:max(0.5));
  // 100 ms
  // maxAttack = 0.1;
  // 2 sec
  // maxRelease = 2;

  // https://www.desmos.com/calculator/qcjwfaaqc5
