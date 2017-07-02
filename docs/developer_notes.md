# ddf-pipeline developer notes

## general

We are happy to accept pull requests that stick to the general coding
style used in the current scripts. (You will notice PEP-8 is not
strictly adhered to.)

## adding a parameter

If you want to add a parameter, edit `utils/parset.py` to give it a
name and description following the format of the other entries in the
parameter description tuple. Parameters are available to the main code
in the options dictionary `o`.
