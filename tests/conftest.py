import os

# The utils/ modules import each other as a package (``from utils.x import ...``)
# only when DDF_LOCAL_DEV=1. Set it before any test imports the package.
os.environ["DDF_LOCAL_DEV"] = "1"
