declare name "lamb";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPLv3";

import("stdfaust.lib");


// a lin needs:
// N times:
// idX,prevSize, midX
// one time:
// offSet, wf, mid

lin2d =
  it.interpolate_linear(
    dy
  , it.interpolate_linear(dx,v0,v1)
  , it.interpolate_linear(dx,v2,v3))
with {
  i0 = rid(int(idX), midX, C)+yOffset;
  i1 = i0+1;
  i2 = i0+sizeX;
  i3 = i1+sizeX;
  dx  = idX-int(idX);
  dy  = idY-int(idY);
  v0 = rdtable(size, wf, rid(i0, mid, C));
  v1 = rdtable(size, wf, rid(i1, mid, C));
  v2 = rdtable(size, wf, rid(i2, mid, C));
  v3 = rdtable(size, wf, rid(i3, mid, C));
};
N = 2;
il(v0,v1,x) = it.interpolate_linear(x,v0,v1);
// inputs: dv,v0,v1
lin(1) = it.interpolate_linear;
// inputs for N=2:
// dy dx v0 v1 dx v2 v3
// inputs for N=3:
// dz dy dx v0 v1 dx v2 v3 dy dx v4 v5 dx v6 v7
lin(N) =
  (_,
   (
     // ((_<:(_,_)),si.bus(prevNrIn*2-2))
     // : (_,ro.crossNM(1,prevNrIn-1),si.bus(prevNrIn-1))
     // :
     ( (si.bus(prevNrIn)<:si.bus(prevNrIn*2)) , si.bus(prevNrIn))
     : (si.bus(prevNrIn),ro.crossnn(prevNrIn) ))
  )
  :(ro.crossNM(1,prevNrIn)
   ,si.bus(prevNrIn*2)
    // , (!,si.bus(prevNrIn-1))
   )
  : il(lin(N-1),lin(N-1),_)
with {
  prevNrIn = inputs(lin(N-1));
};
process =
  // lin(3);
  // tabulateNd(N,0,pwrSine,sizeX,sizeY,rx0,ry0,rx1,ry1,x,y)
  // tabulate2d(0,pwrSine,sizeX,sizeY,rx0,ry0,rx1,ry1,x,y).val(x,y)
  // , pwrSine(x,y);
  tabulateNd(3,0,pwrSineDiv,sizeX,sizeY,sizeY,rx0,ry0,0,rx1,ry1,1,x,y,z)
, pwrSineDiv(x,y,z);

tabulateNd(N,C,expression) =
  calc
  // .readIndex
  // .wf
  // , calc.ri
  .val
with {
  calc =
    environment {
      // total size of the table: s(0) * s(1)  ...  * s(N-2) * s(N-1)
      // N in, 1 out
      size(1) = _;
      size(N) = _*size(N-1);
      // Maximum indexes to access
      // N in, N out
      mids = par(i, N, _-1);
      // Maximum total index to access
      mid = size(N)-1;
      // Prepare the 'float' table read index for one parameter
      idp(midX,r0,r1,x) = (x-r0)/(r1-r0)*midX;
      // Prepare the 'float' table read index for all parameters
      ids =
        ro.interleave(N,4)
        : par(i, N, idp) ;

      // one waveform parameter write value:
      wfp(prevSize,midX,sizeX,r0,r1) =
        r0+
        ((float(
             floor(ba.time%(prevSize*sizeX)/prevSize)
           )*(r1-r0)
         )
         /float(midX))
       ,(prevSize*sizeX);

      // all waveform parameters write values:
      wfps =
        (
          // from sizes
          (bs<:(mids,bs))
          // from r0 r1
        , si.bus(N*2)
        )
        :ro.interleave(N,4)
        : (1,si.bus(4*N))
        : seq(i, N, si.bus(i),wfp, si.bus(4*N-(4*(i+1))))
        : (si.bus(N),!)
      ;

      // Create the table
      wf = wfps:expression;

      // Limit the table read index in [0, mid] if C = 1
      rid(x,mid, 0) = x;
      rid(x,mid, 1) = max(0, min(x, mid));

      // Tabulate an unary 'FX' function on a range [r0, r1]
      val =
        ( (si.bus(3*N)<:si.bus(6*N)), bs )
        :
        ( (bs<:si.bus(2*N)) , si.bus(6*N) )
        // :
        : (size(N),wf,readIndex)
        : rdtable
      ;
      readIndex
      // (sizes,r0s,r1s,xs)
      =
        // from sizes
        (midSizesMidsFromSizes
         // from r0s,r1s,xs
        , si.bus(N*3))
        : (
        // from sizes, mid
        si.bus(N+1)
        // from mids:
      , ( bs<:si.bus(N*2) )
        // from r0s,r1s,xs
      , si.bus(N*3)
      ) :
        // from sizes, mid,mids
        (si.bus(1+2*N)
        ,ids)  // takes (midX,r0,r1,x)
        // output:
        // mid, sizes, mids, ids
        : ri
          // output:
          // mid, total size, read index
        : riPost
      ;
      riPost(mid,size,ri) =
        rid(ri,mid,C);
      midSizesMidsFromSizes = bs<:(mid,bs,mids);
      ri =
        // mid
        _
      , (
        ro.interleave(N,3)
        : (1,0,si.bus(3*N))
        : seq(i, N, riN, si.bus(3*(N-i-1)))
      );

      riN(prevSize,prevID,sizeX,midX,idX) =
        (prevSize*sizeX)
      , ( (prevSize*
           rid(floor(idX),midX,C))
          +prevID) ;

      // shortcut
      bs = si.bus(N);
    };
};

tabulate2d(C,expression,sizeX,sizeY, rx0, ry0, rx1, ry1,x,y) =
  environment {
    size = sizeX*sizeY;
    // Maximum X index to access
    midX = sizeX-1;
    // Maximum Y index to access
    midY = sizeY-1;
    // Maximum total index to access
    mid = size-1;
    // Create the table
    wf = expression(wfX,wfY);
    // Prepare the 'float' table read index for X
    idX = (x-rx0)/(rx1-rx0)*midX;
    // Prepare the 'float' table read index for Y
    idY = ((y-ry0)/(ry1-ry0))*midY;
    // table creation X:
    wfX =
      rx0+float(ba.time%sizeX)*(rx1-rx0)
      /float(midX);
    // table creation Y:
    wfY =
      ry0+
      ((float(ba.time-(ba.time%sizeX))
        /float(sizeX))
       *(ry1-ry0)
      )
      /float(midY)
    ;

    // Limit the table read index in [0, mid] if C = 1
    rid(x,mid, 0) = x;
    rid(x,mid, 1) = max(0, min(x, mid));

    // Tabulate an unary 'FX' function on a range [r0, r1]
    val(x,y) =
      rdtable(size, wf, readIndex);
    readIndex =
      rid(
        rid(int(idX),midX, C)
        +yOffset
      , mid, C);
    yOffset =
      sizeX*rid(floor(idY),midY,C);

    // Tabulate an unary 'FX' function over the range [r0, r1] with linear interpolation
    lin =
      it.interpolate_linear(
        dy
      , it.interpolate_linear(dx,v0,v1)
      , it.interpolate_linear(dx,v2,v3))
    with {
      i0 = rid(int(idX), midX, C)+yOffset;
      i1 = i0+1;
      i2 = i0+sizeX;
      i3 = i1+sizeX;
      dx  = idX-int(idX);
      dy  = idY-int(idY);
      v0 = rdtable(size, wf, rid(i0, mid, C));
      v1 = rdtable(size, wf, rid(i1, mid, C));
      v2 = rdtable(size, wf, rid(i2, mid, C));
      v3 = rdtable(size, wf, rid(i3, mid, C));
    };
    // Tabulate an unary 'FX' function over the range [r0, r1] with cubic interpolation
    cub =
      it.interpolate_cubic(
        dy
      , it.interpolate_cubic(dx,v0,v1,v2,v3)
      , it.interpolate_cubic(dx,v4,v5,v6,v7)
      , it.interpolate_cubic(dx,v8,v9,v10,v11)
      , it.interpolate_cubic(dx,v12,v13,v14,v15)
      )
    with {
      i0  = i4-sizeX;
      i1  = i5-sizeX;
      i2  = i6-sizeX;
      i3  = i7-sizeX;

      i4  = i5-1;
      i5  = rid(int(idX), midX, C)+yOffset;
      i6  = i5+1;
      i7  = i6+1;

      i8  = i4+sizeX;
      i9  = i5+sizeX;
      i10 = i6+sizeX;
      i11 = i7+sizeX;

      i12 = i4+(2*sizeX);
      i13 = i5+(2*sizeX);
      i14 = i6+(2*sizeX);
      i15 = i7+(2*sizeX);

      dx  = idX-int(idX);
      dy  = idY-int(idY);
      v0  = rdtable(size, wf, rid(i0 , mid, C));
      v1  = rdtable(size, wf, rid(i1 , mid, C));
      v2  = rdtable(size, wf, rid(i2 , mid, C));
      v3  = rdtable(size, wf, rid(i3 , mid, C));
      v4  = rdtable(size, wf, rid(i4 , mid, C));
      v5  = rdtable(size, wf, rid(i5 , mid, C));
      v6  = rdtable(size, wf, rid(i6 , mid, C));
      v7  = rdtable(size, wf, rid(i7 , mid, C));
      v8  = rdtable(size, wf, rid(i8 , mid, C));
      v9  = rdtable(size, wf, rid(i9 , mid, C));
      v10 = rdtable(size, wf, rid(i10, mid, C));
      v11 = rdtable(size, wf, rid(i11, mid, C));
      v12 = rdtable(size, wf, rid(i12, mid, C));
      v13 = rdtable(size, wf, rid(i13, mid, C));
      v14 = rdtable(size, wf, rid(i14, mid, C));
      v15 = rdtable(size, wf, rid(i15, mid, C));
    };
  };

sineShaper(x) = (sin((x*.5 + 0.75)*2*ma.PI)+1)*0.5;
pwr(x) = pow(2,x);
pwrSine(x,y)=
  sineShaper(x *(1+(y/ry1))) ;
pwrSineDiv(x,y,z) = pwrSine(x,y)/(1+z);


// x = (float((hslider("x", 0.2, 0.2, 2, 0.01)/2)*midX:floor)*2.0)/midX;
// y = (float((hslider("y", 0.3, 0.3, 3, 0.01)/3)*midY:floor)*3.0)/midY;
// x = ((hslider("x", 0, 0, 2, 0.01)/2)*midX:floor/midX)*2;
xs = hslider("x", rx0, rx0, rx1, 0.01)*midX:floor/midX;
xr = (((((hslider("x", rx0, rx0, rx1, 0.01)
          -rx0
         )
         /(rx1-rx0)
        )*midX:floor/midX)
       *(rx1-rx0)
      )
      +rx0)
     *(rx1-rx0)

;
// x= hslider("x", rx0, rx0, rx1, 1.0/sizeX)*midX:floor/midX;
x= hslider("x", rx0, rx0, rx1, 0.01);
y= hslider("y", ry0, ry0, ry1, 0.01);
z= hslider("z", 0, 0, 1, 0.01);
// idX = (x-rx0)/(rx1-rx0)*midX;
rx0 = 0.1;
rx1 = 1.0;
ry0 = 0.3;
ry1 = 0.7;
// y = hslider("y", , 0, 1, 0.01)*midY:floor/midY;
// y = (float((hslider("y", 0, 0, 1, 0.01)/1.0)*midY:floor)*1.0)/midY;
sizeX = 1<<9;
sizeY = 1<<9;
midX = sizeX-1;
midY = sizeY-1;
// process =
oldProc =
  tabulateNd(N,0,pwrSine)
  // tabulate2d(0,pwrSine,sizeX,sizeY,rx0,ry0,rx1,ry1,x,y).val
  // , tabulate2d(0,pwrSine,sizeX,sizeY,rx0,ry0,rx1,ry1,x,y).lin
  // , tabulate2d(0,pwrSine,sizeX,sizeY,rx0,ry0,rx1,ry1,x,y).cub
  // , pwrSine(x,y)
  // hgroup("",
  // vgroup("[2]test", test)
  // :vgroup("[1]AR",
  // AR
  // ))

  // test
  // ARtest<:
  // PMI_FBFFcompressor_N_chan(strength,thresh,att,rel,knee,prePost,link,FBFF,meter,N)
  // (ARtest:PMI_compression_gain_mono_db(strength,thresh,att,rel,knee,prePost):ba.db2linear)
  // , os.lf_sawpos(1)>0.5
;
mysel(x)=         (checkbox("AR")*
                   (si.onePoleSwitching(hslider("rel simple", 8, 0, 1000, 0.1)*0.001
                                       ,hslider("att simple", 8, 0, 1000, 0.1)*0.001,x))
                  , (1-checkbox("AR"))*AR(x):(!,_,!)):>_,x;


AR = loop~(_,_)
          // :(!,si.bus(4))
          // :(!,_)
with {
  loop(prevRamp,prevGain,x) =
    ramp
  , gain
  , x
  , (x:si.onePoleSwitching(releaseOP,attackOP))
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
      shapeDif(shape,ramp,duration,ma.SR)*fullDif;
    fullDif =dif/(1-warpedSine(shape,ramp));
  };
  shapeDif(shape,phase,duration,sr) =
    shapeDifFormula(shape,phase,duration,48000);
  shapeDifFormula(shape,phase,duration,sr) =
    // warpedSineFormula(shape,phase+(1 / sr / duration))
    // - warpedSineFormula(shape,phase);
    // shapeDif(shape,phase,duration,sr) =
    warpedSine(shape,phase+(1 / sr / duration))
    - warpedSine(shape,phase);

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
      shapeDif(shape,x,duration,ma.SR)
      *(1/(1-warpedSine(shape,x)));
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
  , shapeDif(shape,prevRamp+rampStep,duration',ma.SR)
    * ((dif'/dif)/(1-warpedSine(shape',prevRamp)))
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
  warpedSine(shape,x) =
    // at low number of compares the raw formula is faster than the tabulated version
    // 16 compares: 5 to 6 % CPU
    // when doing contant phase recalculations, we need higher precision in the newramp function
    // cause we get wrong ramp durations (to steep or not steep enough) otherwise
    // 21 compares seems to work well enough in all cases so far
    // at the higher number of compares (21) we get 11-12% CPU for the raw formaula
    // warpedSineFormula(shape,x)
    // the tables do much better
    par(i, nrShapes+1, table(i) * xfadeSelector(shapeSlider,i)):>_
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
    SIZE = 1<<16;
    table(i) = ba.tabulate(0, warpedSineFormula(shapeSliderVal(i)), SIZE, 0, 1, x).lin;
    // 16 compares: 4.5 -6%CPU
    // 21 compares: 7 % CPU
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

  shape =
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
