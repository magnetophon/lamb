## lamb


A look ahead compressor/limiter that's soft as a lamb. 

The secret sauce is all in the attack/release,
you can change both the length and the shape of their curve.

The shapes look like [this](https://www.desmos.com/calculator/iuvx0mrsyi); t in desmos corresponds to the shape parameter in the plugin.

When it is in the middle value, the curve is a part of a pure sine. 
