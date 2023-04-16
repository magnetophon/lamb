declare name "lamb";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPLv3";

import("stdfaust.lib");

process =
  // PMI_FBFFcompressor_N_chan(strength,thresh,att,rel,knee,prePost,link,FBFF,meter,N);
  (ARtest:PMI_compression_gain_mono_db(strength,thresh,att,rel,knee,prePost):ba.db2linear)
, os.lf_sawpos(1)>0.5
  // , (os.lf_sawpos(1)>0.5:moog_AR(att,rel))
;

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
  PMI(PMItime) : ba.linear2db : gain_computer(strength,thresh,knee) : ba.bypass1(prePost,ba.db2linear:si.onePoleSwitching(rel,att):ba.linear2db)  : ba.bypass1((1-prePost),si.onePoleSwitching(rel,att))
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

// my_AR(att,rel) = si.onePoleSwitching(a,r) with {
// a =
// };


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
