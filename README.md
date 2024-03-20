# lamb 🐑

A lookahead compressor/limiter that's soft as a lamb. 

There is also a version with a Rust GUI, latency compensation and proper scaling of the parameters:
https://github.com/magnetophon/lamb-rs

Lamb was made with these goals in mind:
- Be as clean as possible
- Give the user full control over the character with the minimum amount of knobs.

The secret sauce is all in the attack/release:
you can change both the length and the shape of their curve.  
The shapes look like [this](https://www.desmos.com/calculator/cog4ujr7cs); _c0_ in Desmos corresponds to the _shape_ parameter in the plugin.  
When it is at value 0, the curve is a slice of pure sine.  

The ``release hold`` parameter prevents the gain reduction from coming back up if it needs to go down again soon.
You control how soon is soon with ``release hold``.
This adds latency though.


## user preferences

  **ATTENTION** If you want to use the plugin with a samplerate of more than 48k, make sure you change 
  MaxSampleRate at the start of lamb.dsp.  
  There's a couple of other user preferences as well, documented in the dsp file.

🐑
